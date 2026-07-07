#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos > 0L) {
      out[[substr(kv, 1L, pos - 1L)]] <- substr(kv, pos + 1L, nchar(kv))
    } else {
      out[[kv]] <- "TRUE"
    }
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) > 1L) return(x)
  if (is.na(x) || !nzchar(trimws(as.character(x)))) y else x
}

as_chr <- function(x, default = "") as.character((x %||% default)[[1L]])

as_bool <- function(x, default = FALSE) {
  y <- tolower(trimws(as_chr(x, if (default) "TRUE" else "FALSE")))
  if (!nzchar(y)) return(default)
  y %in% c("true", "t", "1", "yes", "y", "on")
}

as_int <- function(x, default = NA_integer_) {
  y <- suppressWarnings(as.integer(as_chr(x, as.character(default))))
  if (!is.finite(y) || is.na(y)) default else y
}

as_num <- function(x, default = NA_real_) {
  y <- suppressWarnings(as.numeric(as_chr(x, as.character(default))))
  if (!is.finite(y) || is.na(y)) default else y
}

as_char_vec <- function(x, default = character()) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return(default)
  vals <- trimws(unlist(strsplit(paste(as.character(x), collapse = ","), ",", fixed = TRUE), use.names = FALSE))
  vals <- vals[nzchar(vals)]
  if (length(vals)) vals else default
}

read_tsv <- function(path) {
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_csv_plain <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, quote = FALSE, row.names = FALSE)
  invisible(path)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

rbind_fill_plain <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && is.data.frame(x), rows)
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

sanitize_label <- function(x, fallback = "warmup") {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", trimws(as.character(x)))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) fallback else x
}

local_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

default_project_root <- function() {
  normalizePath(file.path(local_script_dir(), "../../../../.."), mustWork = FALSE)
}

default_result_root <- function(project_root = default_project_root()) {
  file.path(project_root, "oxygen", "results", "analysis", "parameter_landscape")
}

default_invivo_parameter_table <- function(project_root = default_project_root()) {
  file.path(project_root, "oxygen", "data", "O2_supply_demand", "parameter_table_O2.csv")
}

default_invitro_parameter_table <- function(project_root = default_project_root()) {
  file.path(project_root, "oxygen", "data", "O2_supply_demand", "parameter_table_invitro_buffering.csv")
}

metric_map <- function(path) {
  tab <- read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Summary file is not a metric/value table: ", path, call. = FALSE)
  }
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

metric_value <- function(metrics, keys, default = NA_character_) {
  hit <- keys[keys %in% names(metrics)]
  if (!length(hit)) return(default)
  val <- metrics[[hit[[1L]]]]
  if (length(val) == 0L || is.na(val) || !nzchar(trimws(as.character(val)))) default else val
}

metric_numeric <- function(metrics, keys, seed, label) {
  val <- suppressWarnings(as.numeric(metric_value(metrics, keys, default = NA_character_)))
  if (!is.finite(val) || is.na(val)) {
    warning("Missing/non-finite ", label, " in seed", seed)
    return(NA_real_)
  }
  val
}

seed_from_dir <- function(path) {
  as.integer(sub("^seed", "", basename(path)))
}

list_seed_dirs <- function(input_dir) {
  dirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^seed[0-9]+$", basename(dirs))]
  dirs <- dirs[file.exists(file.path(dirs, "fit_summary.tsv")) & file.exists(file.path(dirs, "best_params.tsv"))]
  dirs[order(seed_from_dir(dirs))]
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path, call. = FALSE)
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter/value columns: ", path, call. = FALSE)
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  stats::setNames(vals, as.character(tab$parameter))
}

invitro_parameter_transform_map <- function() {
  data.frame(
    param_symbol = c(
      "lam_max", "p_misseg", "k_o_mis",
      "buffer_smax", "buffer_beta", "buffer_n_exp",
      "p_wgd", "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
      "O2_crit", "n_O", "p_mis_base", "sigma_growth", "sigma_kary",
      "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
    ),
    param_name = c(
      "log10_lam_max", "log10_p_misseg", "log10_k_o_mis",
      "buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp",
      "log10_p_wgd", "log10_alpha_o2", "gamma_growth", "log10_mu_hp", "gamma_mu",
      "log10_O2_crit", "n_O", "log10_p_mis_base", "log10_sigma_growth", "log10_sigma_kary",
      "init_mean_2N", "log10_init_sd_2N", "init_mean_4N", "log10_init_sd_4N"
    ),
    transform = c(
      "log10", "log10", "log10",
      "identity", "log10", "log10",
      "log10", "log10", "identity", "log10", "identity",
      "log10", "identity", "log10", "log10", "log10",
      "identity", "log10", "identity", "log10"
    ),
    stringsAsFactors = FALSE
  )
}

invivo_parameter_transform_map <- function() {
  data.frame(
    param_symbol = c(
      "lam_max", "p_mis_base", "p_misseg", "k_o_mis",
      "buffer_smax", "buffer_beta", "buffer_n_exp",
      "p_wgd", "o2_S0", "kappa_O", "eta_o2", "rho_2N",
      "alpha_o2", "gamma_growth", "mu_hp", "gamma_mu",
      "O2_crit", "n_O", "tau_O2", "k_clear", "sigma_burden"
    ),
    param_name = c(
      "log10_lam_max", "log10_p_mis_base", "log10_p_misseg", "log10_k_o_mis",
      "buffer_smax", "log10_buffer_beta", "log10_buffer_n_exp",
      "log10_p_wgd", "log10_o2_S0", "log10_kappa_O", "log10_eta_o2", "log10_rho_2N",
      "log10_alpha_o2", "gamma_growth", "log10_mu_hp", "gamma_mu",
      "log10_O2_crit", "n_O", "log10_tau_O2", "log10_k_clear", "log10_sigma_burden"
    ),
    transform = c(
      "log10", "log10", "log10", "log10",
      "identity", "log10", "log10",
      "log10", "log10", "log10", "log10", "log10",
      "log10", "identity", "log10", "identity",
      "log10", "identity", "log10", "log10", "log10"
    ),
    stringsAsFactors = FALSE
  )
}

transform_param_slot <- function(value, transform, param_symbol, slot_label) {
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || is.na(value)) stop("Non-finite ", slot_label, " for parameter ", param_symbol, call. = FALSE)
  if (identical(transform, "log10")) {
    if (value <= 0) stop(param_symbol, " ", slot_label, " must be > 0 for log10 transform.", call. = FALSE)
    return(log10(value))
  }
  value
}

convert_symbol_parameter_table <- function(tab, path, transform_map, label) {
  req <- c("param_symbol", "init_value", "lower_bound", "upper_bound")
  missing <- setdiff(req, names(tab))
  if (length(missing)) stop("parameter_table_input.csv missing columns: ", paste(missing, collapse = ", "), " in ", path, call. = FALSE)
  estimate_col <- if ("use_invitro_fit" %in% names(tab)) "use_invitro_fit" else "estimate"
  if (!estimate_col %in% names(tab)) stop("parameter_table_input.csv must contain use_invitro_fit or estimate in ", path, call. = FALSE)
  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  map <- transform_map
  tab <- merge(map, tab, by = "param_symbol", all.x = FALSE, all.y = FALSE, sort = FALSE)
  tab <- tab[match(map$param_symbol, tab$param_symbol, nomatch = 0L), , drop = FALSE]
  if (!nrow(tab)) stop("No supported ", label, " optimizer parameters found in ", path, call. = FALSE)
  data.frame(
    param_name = as.character(tab$param_name),
    estimate = vapply(tab[[estimate_col]], as_bool, logical(1), default = FALSE),
    init_value = mapply(transform_param_slot, tab$init_value, tab$transform, tab$param_symbol, "init_value"),
    lower_bound = mapply(transform_param_slot, tab$lower_bound, tab$transform, tab$param_symbol, "lower_bound"),
    upper_bound = mapply(transform_param_slot, tab$upper_bound, tab$transform, tab$param_symbol, "upper_bound"),
    param_prototype = as.character(tab$param_symbol),
    stringsAsFactors = FALSE
  )
}

convert_invitro_parameter_table <- function(tab, path) {
  convert_symbol_parameter_table(tab, path, invitro_parameter_transform_map(), "in vitro")
}

convert_invivo_parameter_table <- function(tab, path) {
  convert_symbol_parameter_table(tab, path, invivo_parameter_transform_map(), "in vivo")
}

read_param_table <- function(seed_dir, dataset = "invivo", fallback_path = "") {
  candidates <- c(
    file.path(seed_dir, "parameter_table.csv"),
    file.path(dirname(seed_dir), "parameter_table.csv"),
    file.path(seed_dir, "parameter_table_input.csv"),
    file.path(dirname(seed_dir), "parameter_table_input.csv"),
    file.path(seed_dir, "parameter_table_input_invitro.csv"),
    file.path(dirname(seed_dir), "parameter_table_input_invitro.csv"),
    fallback_path
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Missing parameter table for ", seed_dir, call. = FALSE)
  path <- existing[[1L]]
  tab <- read_csv_plain(path)
  req <- c("param_name", "estimate", "init_value", "lower_bound", "upper_bound", "param_prototype")
  missing <- setdiff(req, names(tab))
  if (length(missing)) {
    if ("param_symbol" %in% names(tab)) {
      if (identical(dataset, "invitro")) return(convert_invitro_parameter_table(tab, path))
      return(convert_invivo_parameter_table(tab, path))
    }
    stop("parameter_table.csv missing columns: ", paste(missing, collapse = ", "), " in ", path, call. = FALSE)
  }
  tab$estimate <- vapply(tab$estimate, as_bool, logical(1), default = FALSE)
  tab$init_value <- suppressWarnings(as.numeric(tab$init_value))
  tab$lower_bound <- suppressWarnings(as.numeric(tab$lower_bound))
  tab$upper_bound <- suppressWarnings(as.numeric(tab$upper_bound))
  tab
}

read_np <- function(metrics, param_table) {
  np <- as_int(metric_value(metrics, c("NP", "NP_used", "NP_requested")), NA_integer_)
  if (!is.finite(np) || is.na(np) || np < 1L) {
    np <- max(10L * sum(param_table$estimate), 1L)
    warning("Could not read NP from fit_summary.tsv; using ", np)
  }
  np
}

sample_uniform_box <- function(n, lower, upper) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) return(matrix(numeric(0), nrow = 0L, ncol = d))
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  span <- as.numeric(upper - lower)
  out <- sweep(u, 2L, span, `*`)
  out <- sweep(out, 2L, as.numeric(lower), `+`)
  colnames(out) <- names(lower)
  out
}

build_de_initialpop <- function(np, lower, upper, init_use, mode = "uniform") {
  np <- as.integer(np)
  if (is.na(np) || np < 1L) stop("NP must be >= 1", call. = FALSE)
  mode <- tolower(trimws(as.character(mode)))
  if (!identical(mode, "uniform")) stop("Only uniform initial population mode is supported here.", call. = FALSE)
  pop <- sample_uniform_box(np, lower, upper)
  init_vec <- as.numeric(init_use[names(lower)])
  names(init_vec) <- names(lower)
  if (any(!is.finite(init_vec))) stop("Initial vector has missing/non-finite values", call. = FALSE)
  pop[1L, ] <- pmin(pmax(init_vec, lower), upper)
  if (np > 1L) pop[2L:np, ] <- sample_uniform_box(np - 1L, lower, upper)
  pop
}

inverse_transform_column <- function(param_name, x) {
  if (startsWith(param_name, "log10_")) return(10^as.numeric(x))
  as.numeric(x)
}

initial_population_natural <- function(seed_dir, seed, metrics, best_names, dataset = "invivo", parameter_table_fallback = "") {
  param_table <- read_param_table(seed_dir, dataset = dataset, fallback_path = parameter_table_fallback)
  opt_tab <- param_table[param_table$estimate, , drop = FALSE]
  if (!nrow(opt_tab)) stop("No estimated parameters in parameter table for ", seed_dir, call. = FALSE)
  param_names <- as.character(opt_tab$param_name)
  init <- stats::setNames(as.numeric(opt_tab$init_value), param_names)
  lower <- stats::setNames(as.numeric(opt_tab$lower_bound), param_names)
  upper <- stats::setNames(as.numeric(opt_tab$upper_bound), param_names)
  np <- read_np(metrics, param_table)

  set.seed(as.integer(seed))
  pop_t <- if (identical(dataset, "invitro")) {
    sample_uniform_box(np, lower, upper)
  } else {
    build_de_initialpop(np = np, lower = lower, upper = upper, init_use = init, mode = "uniform")
  }

  out <- data.frame(matrix(nrow = nrow(pop_t), ncol = 0L))
  for (j in seq_along(param_names)) {
    natural_name <- as.character(opt_tab$param_prototype[[j]])
    out[[natural_name]] <- inverse_transform_column(param_names[[j]], pop_t[, j])
  }
  if ("c_vol_2N_eff_mm3" %in% best_names && "rho_2N" %in% names(out)) {
    out[["c_vol_2N_eff_mm3"]] <- 1 / out[["rho_2N"]]
  }
  best_vals <- read_best_params(seed_dir)
  for (nm in setdiff(best_names, names(out))) {
    out[[nm]] <- if (nm %in% names(best_vals)) as.numeric(best_vals[[nm]]) else NA_real_
  }
  out <- out[, best_names, drop = FALSE]
  out$seed <- as.integer(seed)
  out
}

negative_metric_numeric <- function(metrics, keys, seed, label) {
  val <- metric_numeric(metrics, keys, seed, label)
  if (is.finite(val)) -val else NA_real_
}

best_params_row <- function(seed_dir, seed, best_names, dataset = "invivo") {
  vals <- read_best_params(seed_dir)
  out <- as.data.frame(as.list(stats::setNames(rep(NA_real_, length(best_names)), best_names)))
  for (nm in intersect(best_names, names(vals))) out[[nm]] <- as.numeric(vals[[nm]])
  metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
  out$objective <- metric_numeric(
    metrics,
    c("objective", "objective_total", "optimizer_local_objective", "optimizer_deoptim_objective"),
    seed,
    "objective"
  )
  if (identical(dataset, "invivo")) {
    out$objective_ploidy <- metric_numeric(metrics, "objective_ploidy", seed, "objective_ploidy")
    out$objective_burden <- metric_numeric(metrics, "objective_burden", seed, "objective_burden")
  } else {
    out$objective_growth <- negative_metric_numeric(metrics, "growth_loglik", seed, "objective_growth")
    out$objective_ploidy <- negative_metric_numeric(metrics, "ploidy_loglik", seed, "objective_ploidy")
    out$objective_flow <- negative_metric_numeric(metrics, "flow_loglik", seed, "objective_flow")
  }
  iter_completed <- suppressWarnings(as.integer(metric_value(
    metrics,
    c("optimizer_iter_completed", "deoptim_iter_completed", "deoptim_stop_iteration", "iter_completed"),
    default = NA_character_
  )))
  out$deoptim_iter_completed <- iter_completed
  out$seed <- as.integer(seed)
  out
}

build_seed_parameter_tables <- function(dataset, input_dir, tables_dir, max_seeds = NA_integer_, parameter_table_fallback = "") {
  input_dir <- normalizePath(path.expand(input_dir), mustWork = TRUE)
  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No valid seed directories found under: ", input_dir, call. = FALSE)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
  }
  best_vectors <- lapply(seed_dirs, read_best_params)
  best_names <- names(best_vectors[[1L]])
  for (vals in best_vectors[-1L]) best_names <- c(best_names, setdiff(names(vals), best_names))

  message("Found ", length(seed_dirs), " ", dataset, " seed directories.")
  message("Using ", length(best_names), " ", dataset, " natural parameter columns.")
  best_rows <- vector("list", length(seed_dirs))
  init_rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- seed_from_dir(seed_dir)
    metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
    best_rows[[i]] <- best_params_row(seed_dir, seed, best_names, dataset = dataset)
    init_rows[[i]] <- initial_population_natural(seed_dir, seed, metrics, best_names, dataset = dataset, parameter_table_fallback = parameter_table_fallback)
    if (i %% 25L == 0L || i == length(seed_dirs)) message("Processed ", dataset, " ", i, "/", length(seed_dirs), " seeds.")
  }
  best_df <- do.call(rbind, best_rows)
  init_df <- do.call(rbind, init_rows)
  best_path <- file.path(tables_dir, paste0(dataset, "_best_params_by_seed.csv"))
  init_path <- file.path(tables_dir, paste0(dataset, "_deoptim_initial_population.csv"))
  write_csv(best_df, best_path)
  write_csv(init_df[, c(best_names, "seed"), drop = FALSE], init_path)
  list(best_csv = best_path, initial_csv = init_path, seed_dirs = seed_dirs)
}

seed_space_best_names_path <- function(out_dir, dataset) {
  file.path(out_dir, paste0(dataset, "_best_parameter_names.tsv"))
}

seed_space_tasks_path <- function(out_dir, dataset) {
  file.path(out_dir, paste0(dataset, "_seed_space_tasks.tsv"))
}

seed_space_shard_path <- function(out_dir, dataset, shard_type, seed) {
  file.path(out_dir, dataset, shard_type, paste0("seed", as.integer(seed), ".csv"))
}

seed_space_final_paths <- function(tables_dir, dataset) {
  list(
    best_csv = file.path(tables_dir, paste0(dataset, "_best_params_by_seed.csv")),
    initial_csv = file.path(tables_dir, paste0(dataset, "_deoptim_initial_population.csv"))
  )
}

read_best_names_file <- function(path) {
  tab <- read_tsv(path)
  if (!"parameter" %in% names(tab)) stop("Best-name table is missing parameter column: ", path, call. = FALSE)
  as.character(tab$parameter)
}

build_seed_space_tasks <- function(dataset,
                                   input_dir,
                                   out_dir,
                                   tables_dir,
                                   max_seeds = NA_integer_,
                                   parameter_table_fallback = "") {
  dataset <- match.arg(dataset, c("invivo", "invitro"))
  input_dir <- normalizePath(path.expand(input_dir), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(out_dir), mustWork = FALSE)
  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No valid seed directories found under: ", input_dir, call. = FALSE)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
  }
  best_vectors <- lapply(seed_dirs, read_best_params)
  best_names <- names(best_vectors[[1L]])
  for (vals in best_vectors[-1L]) best_names <- c(best_names, setdiff(names(vals), best_names))

  best_names_path <- seed_space_best_names_path(out_dir, dataset)
  write_tsv(
    data.frame(parameter_index = seq_along(best_names), parameter = best_names, stringsAsFactors = FALSE),
    best_names_path
  )
  tasks <- data.frame(
    task_id = seq_along(seed_dirs),
    dataset = dataset,
    seed = vapply(seed_dirs, seed_from_dir, integer(1L)),
    seed_dir = normalizePath(seed_dirs, mustWork = TRUE),
    best_names_path = normalizePath(best_names_path, mustWork = FALSE),
    parameter_table_fallback = normalizePath(path.expand(parameter_table_fallback), mustWork = FALSE),
    best_shard_path = vapply(seed_dirs, function(seed_dir) {
      seed_space_shard_path(out_dir, dataset, "best", seed_from_dir(seed_dir))
    }, character(1L)),
    initial_shard_path = vapply(seed_dirs, function(seed_dir) {
      seed_space_shard_path(out_dir, dataset, "initial", seed_from_dir(seed_dir))
    }, character(1L)),
    stringsAsFactors = FALSE
  )
  tasks$status <- ifelse(file.exists(tasks$best_shard_path) & file.exists(tasks$initial_shard_path), "complete", "not_started")
  tasks_path <- seed_space_tasks_path(out_dir, dataset)
  write_tsv(tasks, tasks_path)
  final_paths <- seed_space_final_paths(tables_dir, dataset)
  write_tsv(
    data.frame(
      key = c("dataset", "input_dir", "tasks_tsv", "best_names_path", "tables_dir", "best_csv", "initial_csv", "seed_count"),
      value = c(dataset, input_dir, tasks_path, best_names_path, tables_dir, final_paths$best_csv, final_paths$initial_csv, nrow(tasks)),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, paste0(dataset, "_seed_space_settings.tsv"))
  )
  message("Wrote ", nrow(tasks), " ", dataset, " seed-space tasks: ", tasks_path)
  invisible(list(tasks_tsv = tasks_path, best_names_path = best_names_path, final_paths = final_paths, seed_dirs = seed_dirs))
}

task_row_by_id <- function(tasks, task_id, lookup_col = "task_id") {
  if (!lookup_col %in% names(tasks)) stop("Task table is missing lookup column: ", lookup_col, call. = FALSE)
  key <- suppressWarnings(as.integer(tasks[[lookup_col]]))
  hit <- which(key == as.integer(task_id))
  if (length(hit) == 1L) return(tasks[hit, , drop = FALSE])
  if (length(hit) > 1L) stop("Task lookup column is not unique for task id ", task_id, call. = FALSE)
  if (as.integer(task_id) >= 1L && as.integer(task_id) <= nrow(tasks)) return(tasks[as.integer(task_id), , drop = FALSE])
  stop("Task id ", task_id, " was not found in task table.", call. = FALSE)
}

run_seed_space_task <- function(tasks_tsv,
                                task_id,
                                task_lookup_column = "task_id",
                                skip_existing = TRUE) {
  tasks <- read_tsv(tasks_tsv)
  row <- task_row_by_id(tasks, task_id = task_id, lookup_col = task_lookup_column)
  dataset <- match.arg(as.character(row$dataset[[1L]]), c("invivo", "invitro"))
  seed <- as.integer(row$seed[[1L]])
  best_path <- as.character(row$best_shard_path[[1L]])
  initial_path <- as.character(row$initial_shard_path[[1L]])
  if (isTRUE(skip_existing) && file.exists(best_path) && file.exists(initial_path)) {
    message("Skipping complete seed-space shard: ", dataset, " seed", seed)
    return(invisible(list(best_shard = best_path, initial_shard = initial_path, skipped = TRUE)))
  }
  seed_dir <- normalizePath(as.character(row$seed_dir[[1L]]), mustWork = TRUE)
  best_names <- read_best_names_file(as.character(row$best_names_path[[1L]]))
  fallback <- as.character(row$parameter_table_fallback[[1L]])
  metrics <- metric_map(file.path(seed_dir, "fit_summary.tsv"))
  best <- best_params_row(seed_dir, seed, best_names, dataset = dataset)
  initial <- initial_population_natural(
    seed_dir,
    seed,
    metrics,
    best_names,
    dataset = dataset,
    parameter_table_fallback = fallback
  )
  write_csv(best, best_path)
  write_csv(initial[, c(best_names, "seed"), drop = FALSE], initial_path)
  message("Wrote seed-space shards for ", dataset, " seed", seed)
  invisible(list(best_shard = best_path, initial_shard = initial_path, skipped = FALSE))
}

collect_seed_space_tables <- function(dataset,
                                      tasks_tsv,
                                      tables_dir) {
  dataset <- match.arg(dataset, c("invivo", "invitro"))
  tasks <- read_tsv(tasks_tsv)
  tasks <- tasks[as.character(tasks$dataset) == dataset, , drop = FALSE]
  if (!nrow(tasks)) stop("No ", dataset, " tasks found in ", tasks_tsv, call. = FALSE)
  missing_best <- tasks$best_shard_path[!file.exists(tasks$best_shard_path)]
  missing_initial <- tasks$initial_shard_path[!file.exists(tasks$initial_shard_path)]
  if (length(missing_best) || length(missing_initial)) {
    stop(
      "Cannot collect seed-space tables; missing best shards: ", length(missing_best),
      ", missing initial shards: ", length(missing_initial),
      if (length(missing_best)) paste0("\nFirst missing best: ", missing_best[[1L]]) else "",
      if (length(missing_initial)) paste0("\nFirst missing initial: ", missing_initial[[1L]]) else "",
      call. = FALSE
    )
  }
  best_names <- read_best_names_file(as.character(tasks$best_names_path[[1L]]))
  best_rows <- lapply(tasks$best_shard_path, read_csv_plain)
  initial_rows <- lapply(tasks$initial_shard_path, read_csv_plain)
  best_df <- rbind_fill_plain(best_rows)
  initial_df <- rbind_fill_plain(initial_rows)
  best_extra <- if (identical(dataset, "invivo")) {
    c("objective", "objective_ploidy", "objective_burden", "deoptim_iter_completed", "seed")
  } else {
    c("objective", "objective_growth", "objective_ploidy", "objective_flow", "deoptim_iter_completed", "seed")
  }
  best_cols <- c(best_names, intersect(best_extra, names(best_df)), setdiff(names(best_df), c(best_names, best_extra)))
  initial_cols <- c(best_names, "seed")
  missing_initial_cols <- setdiff(initial_cols, names(initial_df))
  if (length(missing_initial_cols)) stop("Initial shards are missing columns: ", paste(missing_initial_cols, collapse = ", "), call. = FALSE)
  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  final_paths <- seed_space_final_paths(tables_dir, dataset)
  write_csv(best_df[, best_cols, drop = FALSE], final_paths$best_csv)
  write_csv(initial_df[, initial_cols, drop = FALSE], final_paths$initial_csv)
  write_tsv(
    data.frame(
      key = c("dataset", "tasks_tsv", "seed_count", "best_csv", "initial_csv"),
      value = c(dataset, normalizePath(tasks_tsv, mustWork = FALSE), nrow(tasks), final_paths$best_csv, final_paths$initial_csv),
      stringsAsFactors = FALSE
    ),
    file.path(tables_dir, paste0(dataset, "_seed_space_collect_summary.tsv"))
  )
  message("Collected ", dataset, " seed-space tables: ", final_paths$best_csv, " and ", final_paths$initial_csv)
  invisible(final_paths)
}

umap_parameter_set <- function(dataset) {
  dataset <- match.arg(dataset, c("invivo", "invitro"))
  if (identical(dataset, "invitro")) {
    return(c(
      "O2_crit", "mu_hp", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
      "buffer_n_exp", "n_O", "alpha_o2", "gamma_growth", "lam_max", "p_mis_base",
      "p_wgd", "gamma_mu"
    ))
  }
  c(
    "lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
    "buffer_n_exp", "p_wgd", "o2_S0", "kappa_O", "eta_o2", "alpha_o2",
    "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O", "k_clear"
  )
}

umap_log10_parameter_set <- function(dataset) {
  intersect(
    c(
      "lam_max", "p_mis_base", "p_misseg", "k_o_mis", "buffer_beta", "buffer_n_exp",
      "p_wgd", "o2_S0", "kappa_O", "eta_o2", "alpha_o2", "mu_hp", "O2_crit", "k_clear"
    ),
    umap_parameter_set(dataset)
  )
}

pooled_umap_parameter_set <- function() {
  intersect(umap_parameter_set("invivo"), umap_parameter_set("invitro"))
}

pooled_umap_log10_parameter_set <- function() {
  intersect(pooled_umap_parameter_set(), union(umap_log10_parameter_set("invivo"), umap_log10_parameter_set("invitro")))
}

transform_umap_features <- function(df, params, log10_params) {
  missing <- setdiff(params, names(df))
  if (length(missing)) stop("Input table is missing UMAP parameter columns: ", paste(missing, collapse = ", "), call. = FALSE)
  out <- as.data.frame(lapply(params, function(param) {
    vals <- suppressWarnings(as.numeric(df[[param]]))
    if (param %in% log10_params) {
      if (any(!is.finite(vals) | vals <= 0)) stop("Cannot log10-transform non-positive/non-finite values for parameter: ", param, call. = FALSE)
      vals <- log10(vals)
    }
    vals
  }), check.names = FALSE)
  names(out) <- params
  if (any(!is.finite(as.matrix(out)))) stop("UMAP feature matrix contains non-finite values.", call. = FALSE)
  out
}

drop_parameter_table_initial_rows <- function(initial_df) {
  row_in_seed <- ave(seq_len(nrow(initial_df)), initial_df$seed, FUN = seq_along)
  initial_df[row_in_seed != 1L, , drop = FALSE]
}

prepare_pooled_tables <- function(invivo_best_csv,
                                  invivo_initial_csv,
                                  invitro_best_csv,
                                  invitro_initial_csv,
                                  drop_invivo_parameter_table_initial = TRUE,
                                  drop_invitro_parameter_table_initial = TRUE) {
  params <- pooled_umap_parameter_set()
  log10_params <- pooled_umap_log10_parameter_set()
  read_best <- function(path, dataset) {
    df <- read_csv_plain(path)
    if (!"seed" %in% names(df)) stop("Best-parameter CSV must contain a seed column: ", path, call. = FALSE)
    missing_params <- setdiff(params, names(df))
    if (length(missing_params)) stop(dataset, " best table is missing pooled UMAP parameters: ", paste(missing_params, collapse = ", "), call. = FALSE)
    if (!"objective" %in% names(df)) stop(dataset, " best table is missing objective: ", path, call. = FALSE)
    df$dataset <- dataset
    df$point_type <- "best"
    df$source_group <- paste0(dataset, "_best")
    df$seed <- as.integer(df$seed)
    df
  }
  read_initial <- function(path, dataset, drop_parameter_table_initial) {
    df <- read_csv_plain(path)
    if (!"seed" %in% names(df)) stop("Initial population CSV must contain a seed column: ", path, call. = FALSE)
    missing_params <- setdiff(params, names(df))
    if (length(missing_params)) stop(dataset, " initial table is missing pooled UMAP parameters: ", paste(missing_params, collapse = ", "), call. = FALSE)
    df$source_row_id <- seq_len(nrow(df))
    if (isTRUE(drop_parameter_table_initial)) df <- drop_parameter_table_initial_rows(df)
    df$dataset <- dataset
    df$point_type <- "initial"
    df$source_group <- paste0(dataset, "_initial")
    df$seed <- as.integer(df$seed)
    df$objective <- NA_real_
    df
  }
  invivo_best <- read_best(invivo_best_csv, "invivo")
  invitro_best <- read_best(invitro_best_csv, "invitro")
  invivo_initial <- read_initial(invivo_initial_csv, "invivo", drop_invivo_parameter_table_initial)
  invitro_initial <- read_initial(invitro_initial_csv, "invitro", drop_invitro_parameter_table_initial)
  best_df <- rbind_fill_plain(list(invivo_best, invitro_best))
  initial_df <- rbind_fill_plain(list(invivo_initial, invitro_initial))
  best_df$objective <- suppressWarnings(as.numeric(best_df$objective))
  initial_df$objective <- suppressWarnings(as.numeric(initial_df$objective))
  list(
    params = params,
    log10_params = log10_params,
    initial_df = initial_df,
    best_df = best_df,
    initial_features = transform_umap_features(initial_df, params, log10_params),
    best_features = transform_umap_features(best_df, params, log10_params)
  )
}

standardize_features <- function(x) {
  scaled <- scale(as.matrix(x))
  zero_sd <- !is.finite(attr(scaled, "scaled:scale")) | attr(scaled, "scaled:scale") == 0
  if (any(zero_sd)) stop("Feature columns have zero/non-finite SD after z-score: ", paste(names(x)[zero_sd], collapse = ", "), call. = FALSE)
  scaled
}

normalize_reduction <- function(reduction) {
  reduction <- tolower(trimws(as.character(reduction)))
  aliases <- c(umaps = "umap", tsnes = "tsne", pcas = "pca")
  if (reduction %in% names(aliases)) reduction <- aliases[[reduction]]
  if (!reduction %in% c("umap", "tsne")) stop("reduction must be one of: umap, tsne", call. = FALSE)
  reduction
}

reduction_file_suffix <- function(reduction) normalize_reduction(reduction)

reduction_dir_name <- function(reduction) {
  switch(normalize_reduction(reduction), umap = "UMAPs", tsne = "TSNEs")
}

reduction_coordinate_names <- function(reduction) {
  switch(normalize_reduction(reduction), umap = c("UMAP1", "UMAP2"), tsne = c("tSNE1", "tSNE2"))
}

reduction_axis_labels <- function(reduction) {
  switch(normalize_reduction(reduction), umap = c("UMAP 1", "UMAP 2"), tsne = c("t-SNE 1", "t-SNE 2"))
}

run_umap_embedding <- function(feature_mat, label, umap_seed, n_neighbors, min_dist, n_threads) {
  if (!requireNamespace("uwot", quietly = TRUE)) stop("Required R package is not installed: uwot", call. = FALSE)
  message("Running ", label, " UMAP with ", nrow(feature_mat), " points and ", ncol(feature_mat), " dimensions.")
  set.seed(as.integer(umap_seed))
  emb <- uwot::umap(
    feature_mat,
    n_neighbors = min(as.integer(n_neighbors), nrow(feature_mat) - 1L),
    min_dist = as.numeric(min_dist),
    metric = "euclidean",
    n_threads = as.integer(n_threads),
    n_sgd_threads = 1L,
    ret_model = FALSE,
    verbose = TRUE
  )
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb
}

run_tsne_embedding <- function(feature_mat, label, tsne_seed, perplexity, theta, max_iter) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) stop("Required R package is not installed: Rtsne", call. = FALSE)
  mat <- as.matrix(feature_mat)
  storage.mode(mat) <- "double"
  if (nrow(mat) < 4L) stop("t-SNE requires at least 4 rows for ", label, call. = FALSE)
  perplexity <- min(as_num(perplexity, 30), max(1, floor((nrow(mat) - 1L) / 3L)))
  message("Running ", label, " t-SNE with ", nrow(mat), " points, ", ncol(mat), " dimensions, perplexity ", perplexity, ".")
  set.seed(as.integer(tsne_seed))
  emb <- Rtsne::Rtsne(
    mat,
    dims = 2L,
    perplexity = perplexity,
    theta = as_num(theta, 0.5),
    max_iter = as_int(max_iter, 1000L),
    pca = FALSE,
    check_duplicates = FALSE,
    verbose = TRUE
  )$Y
  colnames(emb) <- c("tSNE1", "tSNE2")
  emb
}

run_reduction_embedding <- function(feature_mat,
                                    reduction,
                                    label,
                                    umap_seed = 123L,
                                    n_neighbors = 80L,
                                    min_dist = 0.1,
                                    n_threads = 1L,
                                    tsne_seed = 123L,
                                    tsne_perplexity = 30,
                                    tsne_theta = 0.5,
                                    tsne_max_iter = 1000L) {
  reduction <- normalize_reduction(reduction)
  if (identical(reduction, "umap")) {
    return(run_umap_embedding(feature_mat, label, umap_seed, n_neighbors, min_dist, n_threads))
  }
  run_tsne_embedding(feature_mat, label, tsne_seed, tsne_perplexity, tsne_theta, tsne_max_iter)
}

build_pooled_plot_data <- function(emb, initial_df, best_df, reduction_label, coord_names) {
  meta <- rbind_fill_plain(list(initial_df, best_df))
  out <- data.frame(
    dataset = as.character(meta$dataset),
    point_type = as.character(meta$point_type),
    source_group = as.character(meta$source_group),
    seed = as.integer(meta$seed),
    objective = suppressWarnings(as.numeric(meta$objective)),
    reduction = reduction_label,
    stringsAsFactors = FALSE
  )
  out[[coord_names[[1L]]]] <- emb[, 1L]
  out[[coord_names[[2L]]]] <- emb[, 2L]
  out[, c(coord_names, setdiff(names(out), coord_names)), drop = FALSE]
}

write_reduction_coordinate_table <- function(reduction, emb, pooled, analysis_root) {
  reduction <- normalize_reduction(reduction)
  coord_names <- reduction_coordinate_names(reduction)
  suffix <- reduction_file_suffix(reduction)
  table_dir <- file.path(analysis_root, "pooled_invivo_invitro", reduction_dir_name(reduction), "Tables", "Full")
  stem <- paste0("pooled_invivo_invitro_initial_vs_best_", suffix)
  plot_data <- build_pooled_plot_data(
    emb = emb,
    initial_df = pooled$initial_df,
    best_df = pooled$best_df,
    reduction_label = paste0("pooled_full_", suffix),
    coord_names = coord_names
  )
  coordinate_path <- file.path(table_dir, paste0(stem, "_coordinates.csv"))
  write_csv(plot_data, coordinate_path)
  coordinate_path
}

relabel_clusters_by_embedding <- function(cluster_num, plot_data, coord_names) {
  centers <- stats::aggregate(
    cbind(.embedding_x, .embedding_y) ~ cluster,
    data = data.frame(
      cluster = as.integer(cluster_num),
      .embedding_x = suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]])),
      .embedding_y = suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
    ),
    FUN = stats::median
  )
  centers <- centers[order(centers$.embedding_x, centers$.embedding_y), , drop = FALSE]
  relabel_map <- stats::setNames(seq_len(nrow(centers)), as.character(centers$cluster))
  as.integer(relabel_map[as.character(cluster_num)])
}

single_cluster_assignment <- function(n, cluster_source, sample_n) {
  list(
    cluster_num = rep(1L, n),
    cluster_id = rep("C01", n),
    summary = data.frame(
      cluster_source = cluster_source,
      k = 1L,
      average_silhouette = NA_real_,
      sample_n = as.integer(sample_n),
      selected = TRUE,
      stringsAsFactors = FALSE
    )
  )
}

auto_silhouette_kmeans <- function(basis_mat,
                                   plot_data,
                                   cluster_source,
                                   coord_names,
                                   seed = 123L,
                                   k_min = 2L,
                                   k_max = 8L,
                                   silhouette_sample_n = 5000L,
                                   nstart = 10L,
                                   iter.max = 100L) {
  if (!requireNamespace("cluster", quietly = TRUE)) stop("Required R package is not installed: cluster", call. = FALSE)
  mat <- as.matrix(basis_mat)
  storage.mode(mat) <- "double"
  if (nrow(mat) != nrow(plot_data)) stop("Cluster matrix row count mismatch for ", cluster_source, call. = FALSE)
  if (!nrow(mat) || !ncol(mat)) stop("Cluster feature matrix is empty for ", cluster_source, call. = FALSE)
  if (any(!is.finite(mat))) stop("Cluster feature matrix contains non-finite values for ", cluster_source, call. = FALSE)

  n <- nrow(mat)
  silhouette_sample_n <- as.integer(silhouette_sample_n)
  if (!is.finite(silhouette_sample_n) || is.na(silhouette_sample_n) || silhouette_sample_n < 2L) {
    silhouette_sample_n <- min(5000L, n)
  }
  sample_n <- min(n, silhouette_sample_n)
  if (n < 3L || sample_n < 3L) return(single_cluster_assignment(n, cluster_source, sample_n))
  k_min <- max(2L, as.integer(k_min))
  k_max <- min(as.integer(k_max), sample_n - 1L)
  set.seed(as.integer(seed))
  sample_idx <- if (n > sample_n) sort(sample.int(n, sample_n)) else seq_len(n)
  sample_mat <- mat[sample_idx, , drop = FALSE]
  unique_sample_n <- nrow(unique(sample_mat))
  k_max <- min(k_max, unique_sample_n - 1L)
  if (!is.finite(k_max) || is.na(k_max) || k_max < k_min) return(single_cluster_assignment(n, cluster_source, sample_n))

  k_values <- seq.int(k_min, k_max)
  sample_dist <- stats::dist(sample_mat)
  scores <- rep(NA_real_, length(k_values))
  centers_by_k <- vector("list", length(k_values))
  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    set.seed(as.integer(seed) + k)
    km <- try(stats::kmeans(sample_mat, centers = k, nstart = as.integer(nstart), iter.max = as.integer(iter.max)), silent = TRUE)
    if (inherits(km, "try-error") || length(unique(km$cluster)) < 2L) next
    sil <- try(cluster::silhouette(km$cluster, sample_dist), silent = TRUE)
    if (inherits(sil, "try-error")) next
    scores[[i]] <- mean(sil[, "sil_width"])
    centers_by_k[[i]] <- km$centers
  }
  summary <- data.frame(
    cluster_source = cluster_source,
    k = k_values,
    average_silhouette = scores,
    sample_n = sample_n,
    selected = FALSE,
    stringsAsFactors = FALSE
  )
  valid <- is.finite(summary$average_silhouette)
  if (!any(valid)) return(single_cluster_assignment(n, cluster_source, sample_n))
  selected_i <- which(valid)[which.max(summary$average_silhouette[valid])]
  selected_k <- summary$k[[selected_i]]
  summary$selected[[selected_i]] <- TRUE
  final_km <- try(stats::kmeans(mat, centers = centers_by_k[[selected_i]], iter.max = as.integer(iter.max), algorithm = "Lloyd"), silent = TRUE)
  if (inherits(final_km, "try-error")) {
    set.seed(as.integer(seed) + selected_k + 1000L)
    final_km <- stats::kmeans(mat, centers = selected_k, nstart = as.integer(nstart), iter.max = as.integer(iter.max))
  }
  cluster_num <- relabel_clusters_by_embedding(final_km$cluster, plot_data, coord_names)
  list(cluster_num = cluster_num, cluster_id = sprintf("C%02d", cluster_num), summary = summary)
}

cluster_dataset_specs <- function() {
  data.frame(
    dataset = c("invivo", "invitro"),
    dataset_label = c("in vivo", "in vitro"),
    output_token = c("invivo", "invitro"),
    cluster_prefix = c("vi", "vt"),
    stringsAsFactors = FALSE
  )
}

summarize_best_clusters <- function(clustered_best, coord_names) {
  rows <- lapply(sort(unique(clustered_best$cluster_id)), function(cluster_id) {
    d <- clustered_best[clustered_best$cluster_id == cluster_id, , drop = FALSE]
    base <- data.frame(
      dataset = unique(d$dataset)[[1L]],
      dataset_label = unique(d$dataset_label)[[1L]],
      cluster_prefix = unique(d$cluster_prefix)[[1L]],
      cluster_id = cluster_id,
      cluster_base_id = unique(d$cluster_base_id)[[1L]],
      cluster_num = unique(d$cluster_num)[[1L]],
      n_seeds = nrow(d),
      seed_min = min(d$seed, na.rm = TRUE),
      seed_max = max(d$seed, na.rm = TRUE),
      objective_mean = mean(d$objective, na.rm = TRUE),
      objective_median = stats::median(d$objective, na.rm = TRUE),
      objective_min = min(d$objective, na.rm = TRUE),
      objective_max = max(d$objective, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    coord_values <- c(
      stats::median(d[[coord_names[[1L]]]], na.rm = TRUE),
      stats::median(d[[coord_names[[2L]]]], na.rm = TRUE),
      min(d[[coord_names[[1L]]]], na.rm = TRUE),
      max(d[[coord_names[[1L]]]], na.rm = TRUE),
      min(d[[coord_names[[2L]]]], na.rm = TRUE),
      max(d[[coord_names[[2L]]]], na.rm = TRUE)
    )
    coord_out <- as.data.frame(
      as.list(stats::setNames(
        coord_values,
        c(
          paste0(coord_names, "_median"),
          paste0(coord_names[[1L]], c("_min", "_max")),
          paste0(coord_names[[2L]], c("_min", "_max"))
        )
      )),
      check.names = FALSE
    )
    cbind(base, coord_out)
  })
  out <- do.call(rbind, rows)
  out[order(out$cluster_num), , drop = FALSE]
}

best_seed_group_table <- function(clustered_best, reduction, coord_names) {
  out <- data.frame(
    method = reduction,
    reduction = reduction,
    dataset = as.character(clustered_best$dataset),
    dataset_label = as.character(clustered_best$dataset_label),
    cluster_prefix = as.character(clustered_best$cluster_prefix),
    cluster_id = as.character(clustered_best$cluster_id),
    cluster_base_id = as.character(clustered_best$cluster_base_id),
    cluster_num = as.integer(clustered_best$cluster_num),
    seed = as.integer(clustered_best$seed),
    objective = suppressWarnings(as.numeric(clustered_best$objective)),
    coordinate_1_name = coord_names[[1L]],
    coordinate_2_name = coord_names[[2L]],
    stringsAsFactors = FALSE
  )
  out$coordinate_1 <- suppressWarnings(as.numeric(clustered_best[[coord_names[[1L]]]]))
  out$coordinate_2 <- suppressWarnings(as.numeric(clustered_best[[coord_names[[2L]]]]))
  out <- out[order(out$method, out$dataset, out$cluster_num, -out$objective, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

cluster_best_dataset <- function(plot_data,
                                 dataset,
                                 dataset_label,
                                 cluster_prefix,
                                 coord_names,
                                 cluster_seed,
                                 cluster_k_min,
                                 cluster_k_max,
                                 silhouette_sample_n) {
  best_idx <- which(plot_data$dataset == dataset & plot_data$point_type == "best")
  if (!length(best_idx)) stop("No rows matched dataset == '", dataset, "' and point_type == 'best'.", call. = FALSE)
  best <- plot_data[best_idx, , drop = FALSE]
  cluster_source <- paste0(dataset, "_best_", paste(coord_names, collapse = "_"))
  assignment <- auto_silhouette_kmeans(
    basis_mat = as.matrix(best[, coord_names, drop = FALSE]),
    plot_data = best,
    cluster_source = cluster_source,
    coord_names = coord_names,
    seed = cluster_seed,
    k_min = cluster_k_min,
    k_max = cluster_k_max,
    silhouette_sample_n = silhouette_sample_n
  )
  selected_summary <- assignment$summary[assignment$summary$selected, , drop = FALSE]
  best$dataset_label <- dataset_label
  best$cluster_scope <- paste0(dataset, "_best")
  best$cluster_source <- cluster_source
  best$cluster_prefix <- cluster_prefix
  best$cluster_base_id <- assignment$cluster_id
  best$cluster_id <- paste0(cluster_prefix, "_", assignment$cluster_id)
  best$cluster_num <- assignment$cluster_num
  best$cluster_k <- length(unique(assignment$cluster_num))
  best$cluster_silhouette_avg <- selected_summary$average_silhouette[[1L]]
  best$cluster_silhouette_sample_n <- selected_summary$sample_n[[1L]]
  sil <- assignment$summary
  sil$dataset <- dataset
  sil$dataset_label <- dataset_label
  sil$cluster_prefix <- cluster_prefix
  sil$cluster_source <- cluster_source
  summary <- summarize_best_clusters(best, coord_names = coord_names)
  summary$cluster_source <- cluster_source
  summary$selected_k <- selected_summary$k[[1L]]
  summary$selected_average_silhouette <- selected_summary$average_silhouette[[1L]]
  list(best_idx = best_idx, best = best, silhouette = sil, summary = summary, selected_summary = selected_summary)
}

seed_key <- function(dataset, seed) paste(as.character(dataset), as.integer(seed), sep = "::")

load_pooled_best_feature_data <- function(invivo_best_csv, invitro_best_csv) {
  params <- pooled_umap_parameter_set()
  log10_params <- pooled_umap_log10_parameter_set()
  read_best <- function(path, dataset) {
    df <- read_csv_plain(path)
    required <- c("seed", params)
    missing <- setdiff(required, names(df))
    if (length(missing)) stop("Best-parameter table for ", dataset, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    df$dataset <- dataset
    df$seed <- suppressWarnings(as.integer(df$seed))
    df
  }
  best_df <- rbind_fill_plain(list(read_best(invivo_best_csv, "invivo"), read_best(invitro_best_csv, "invitro")))
  feature_df <- transform_umap_features(best_df, params, log10_params)
  feature_df <- cbind(data.frame(dataset = as.character(best_df$dataset), seed = as.integer(best_df$seed), stringsAsFactors = FALSE), feature_df)
  feature_df$.seed_key <- seed_key(feature_df$dataset, feature_df$seed)
  list(params = params, log10_params = log10_params, feature_df = feature_df)
}

align_best_feature_rows <- function(clustered_best, feature_data) {
  keys <- seed_key(clustered_best$dataset, clustered_best$seed)
  idx <- match(keys, feature_data$feature_df$.seed_key)
  if (any(is.na(idx))) stop("Missing raw feature rows for clustered best seeds: ", paste(head(keys[is.na(idx)], 20), collapse = ", "), call. = FALSE)
  feature_data$feature_df[idx, feature_data$params, drop = FALSE]
}

zscore_primary_features <- function(raw_features, params, log10_params, method, dataset, primary_cluster_id) {
  mat <- as.matrix(raw_features[, params, drop = FALSE])
  storage.mode(mat) <- "double"
  centers <- colMeans(mat, na.rm = TRUE)
  scales <- apply(mat, 2L, stats::sd, na.rm = TRUE)
  used <- is.finite(centers) & is.finite(scales) & scales > 0
  metadata <- data.frame(
    method = method,
    reduction = method,
    dataset = dataset,
    primary_cluster_id = primary_cluster_id,
    feature = params,
    log10_transformed = params %in% log10_params,
    center = as.numeric(centers),
    scale = as.numeric(scales),
    used = as.logical(used),
    drop_reason = ifelse(used, "", "zero_or_nonfinite_sd"),
    stringsAsFactors = FALSE
  )
  if (!any(used)) return(list(mat = matrix(nrow = nrow(mat), ncol = 0L), metadata = metadata))
  scaled <- scale(mat[, used, drop = FALSE], center = centers[used], scale = scales[used])
  colnames(scaled) <- params[used]
  list(mat = scaled, metadata = metadata)
}

subcluster_seed_group_table <- function(subclustered_best, coord_names) {
  out <- data.frame(
    method = as.character(subclustered_best$method),
    reduction = as.character(subclustered_best$reduction),
    dataset = as.character(subclustered_best$dataset),
    dataset_label = as.character(subclustered_best$dataset_label),
    primary_cluster_id = as.character(subclustered_best$primary_cluster_id),
    primary_cluster_base_id = as.character(subclustered_best$primary_cluster_base_id),
    primary_cluster_num = as.integer(subclustered_best$primary_cluster_num),
    subcluster_id = as.character(subclustered_best$subcluster_id),
    subcluster_base_id = as.character(subclustered_best$subcluster_base_id),
    subcluster_num = as.integer(subclustered_best$subcluster_num),
    seed = as.integer(subclustered_best$seed),
    objective = suppressWarnings(as.numeric(subclustered_best$objective)),
    coordinate_1_name = coord_names[[1L]],
    coordinate_2_name = coord_names[[2L]],
    stringsAsFactors = FALSE
  )
  out$coordinate_1 <- suppressWarnings(as.numeric(subclustered_best[[coord_names[[1L]]]]))
  out$coordinate_2 <- suppressWarnings(as.numeric(subclustered_best[[coord_names[[2L]]]]))
  out <- out[order(out$method, out$dataset, out$primary_cluster_num, out$subcluster_num, -out$objective, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

summarize_best_subclusters <- function(subclustered_best, coord_names) {
  rows <- lapply(sort(unique(subclustered_best$subcluster_id)), function(subcluster_id) {
    d <- subclustered_best[subclustered_best$subcluster_id == subcluster_id, , drop = FALSE]
    min_i <- which.min(d$objective)
    base <- data.frame(
      method = unique(d$method)[[1L]],
      reduction = unique(d$reduction)[[1L]],
      dataset = unique(d$dataset)[[1L]],
      dataset_label = unique(d$dataset_label)[[1L]],
      primary_cluster_id = unique(d$primary_cluster_id)[[1L]],
      primary_cluster_base_id = unique(d$primary_cluster_base_id)[[1L]],
      primary_cluster_num = unique(d$primary_cluster_num)[[1L]],
      subcluster_id = subcluster_id,
      subcluster_base_id = unique(d$subcluster_base_id)[[1L]],
      subcluster_num = unique(d$subcluster_num)[[1L]],
      n_seeds = nrow(d),
      objective_min_seed = d$seed[[min_i]],
      objective_mean = mean(d$objective, na.rm = TRUE),
      objective_median = stats::median(d$objective, na.rm = TRUE),
      objective_min = min(d$objective, na.rm = TRUE),
      objective_max = max(d$objective, na.rm = TRUE),
      subcluster_k = unique(d$subcluster_k)[[1L]],
      subcluster_silhouette_avg = unique(d$subcluster_silhouette_avg)[[1L]],
      subcluster_silhouette_sample_n = unique(d$subcluster_silhouette_sample_n)[[1L]],
      stringsAsFactors = FALSE
    )
    coord_values <- c(
      stats::median(d[[coord_names[[1L]]]], na.rm = TRUE),
      stats::median(d[[coord_names[[2L]]]], na.rm = TRUE),
      min(d[[coord_names[[1L]]]], na.rm = TRUE),
      max(d[[coord_names[[1L]]]], na.rm = TRUE),
      min(d[[coord_names[[2L]]]], na.rm = TRUE),
      max(d[[coord_names[[2L]]]], na.rm = TRUE)
    )
    coord_out <- as.data.frame(
      as.list(stats::setNames(
        coord_values,
        c(
          paste0(coord_names, "_median"),
          paste0(coord_names[[1L]], c("_min", "_max")),
          paste0(coord_names[[2L]], c("_min", "_max"))
        )
      )),
      check.names = FALSE
    )
    cbind(base, coord_out)
  })
  out <- do.call(rbind, rows)
  out[order(out$method, out$dataset, out$primary_cluster_num, out$subcluster_num), , drop = FALSE]
}

subcluster_zscore_feature_table <- function(subclustered_best, zscore_features_by_primary, coord_names) {
  id_cols <- subcluster_seed_group_table(subclustered_best, coord_names = coord_names)
  feature_cols <- unique(unlist(lapply(zscore_features_by_primary, colnames), use.names = FALSE))
  if (!length(feature_cols)) return(id_cols)
  feature_out <- as.data.frame(matrix(NA_real_, nrow = nrow(subclustered_best), ncol = length(feature_cols)))
  names(feature_out) <- paste0("z_", feature_cols)
  cursor <- 1L
  for (feature_mat in zscore_features_by_primary) {
    n <- nrow(feature_mat)
    if (n > 0L && ncol(feature_mat) > 0L) {
      feature_out[cursor:(cursor + n - 1L), paste0("z_", colnames(feature_mat))] <- as.data.frame(feature_mat)
    }
    cursor <- cursor + n
  }
  cbind(id_cols, feature_out)
}

analyze_best_subclusters <- function(clustered_best,
                                     feature_data,
                                     reduction,
                                     coord_names,
                                     subcluster_seed,
                                     subcluster_k_min,
                                     subcluster_k_max,
                                     subcluster_min_n,
                                     silhouette_sample_n) {
  feature_rows <- align_best_feature_rows(clustered_best, feature_data)
  clustered_best$.row_id <- seq_len(nrow(clustered_best))
  split_ids <- split(clustered_best$.row_id, interaction(clustered_best$dataset, clustered_best$cluster_id, drop = TRUE, lex.order = TRUE))
  subclustered <- list()
  silhouettes <- list()
  feature_meta <- list()
  zscore_features <- list()
  for (i in seq_along(split_ids)) {
    idx <- split_ids[[i]]
    d <- clustered_best[idx, , drop = FALSE]
    dataset <- unique(as.character(d$dataset))[[1L]]
    primary_cluster_id <- unique(as.character(d$cluster_id))[[1L]]
    raw <- feature_rows[idx, , drop = FALSE]
    z <- zscore_primary_features(raw, feature_data$params, feature_data$log10_params, reduction, dataset, primary_cluster_id)
    source <- paste0(reduction, "_", primary_cluster_id, "_best_raw_features_zscore")
    if (nrow(d) < as.integer(subcluster_min_n) || ncol(z$mat) == 0L) {
      assignment <- single_cluster_assignment(nrow(d), source, nrow(d))
    } else {
      assignment <- auto_silhouette_kmeans(
        basis_mat = z$mat,
        plot_data = d,
        cluster_source = source,
        coord_names = coord_names,
        seed = as.integer(subcluster_seed) + i - 1L,
        k_min = subcluster_k_min,
        k_max = min(as.integer(subcluster_k_max), nrow(d) - 1L),
        silhouette_sample_n = silhouette_sample_n
      )
    }
    selected_summary <- assignment$summary[assignment$summary$selected, , drop = FALSE]
    d$method <- reduction
    d$reduction <- reduction
    d$primary_cluster_id <- as.character(d$cluster_id)
    d$primary_cluster_base_id <- as.character(d$cluster_base_id)
    d$primary_cluster_num <- as.integer(d$cluster_num)
    d$subcluster_source <- source
    d$subcluster_base_id <- sprintf("Sc%02d", assignment$cluster_num)
    d$subcluster_id <- paste0(d$primary_cluster_id, "_", d$subcluster_base_id)
    d$subcluster_num <- as.integer(assignment$cluster_num)
    d$subcluster_k <- length(unique(assignment$cluster_num))
    d$subcluster_silhouette_avg <- selected_summary$average_silhouette[[1L]]
    d$subcluster_silhouette_sample_n <- selected_summary$sample_n[[1L]]
    d$.row_id <- NULL
    sil <- assignment$summary
    sil$method <- reduction
    sil$reduction <- reduction
    sil$dataset <- dataset
    sil$primary_cluster_id <- primary_cluster_id
    sil$primary_cluster_n <- nrow(d)
    silhouettes[[i]] <- sil[, c("method", "reduction", "dataset", "primary_cluster_id", "primary_cluster_n", setdiff(names(sil), c("method", "reduction", "dataset", "primary_cluster_id", "primary_cluster_n"))), drop = FALSE]
    subclustered[[i]] <- d
    feature_meta[[i]] <- z$metadata
    zscore_features[[i]] <- z$mat
  }
  subclustered_best <- do.call(rbind, subclustered)
  row.names(subclustered_best) <- NULL
  zscore_table <- subcluster_zscore_feature_table(subclustered_best, zscore_features, coord_names = coord_names)
  subclustered_best <- subclustered_best[order(subclustered_best$dataset, subclustered_best$primary_cluster_num, subclustered_best$subcluster_num, subclustered_best$seed), , drop = FALSE]
  row.names(subclustered_best) <- NULL
  zscore_table <- zscore_table[order(zscore_table$dataset, zscore_table$primary_cluster_num, zscore_table$subcluster_num, zscore_table$seed), , drop = FALSE]
  row.names(zscore_table) <- NULL
  list(
    best = subclustered_best,
    seed_groups = subcluster_seed_group_table(subclustered_best, coord_names = coord_names),
    summary = summarize_best_subclusters(subclustered_best, coord_names = coord_names),
    silhouette = do.call(rbind, silhouettes),
    feature_metadata = do.call(rbind, feature_meta),
    zscore_features = zscore_table
  )
}

require_plotting <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required R package is not installed for cluster figures: ggplot2", call. = FALSE)
  }
}

require_pooled_objective_plotting <- function() {
  require_plotting()
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Required R package is not installed for pooled two-scale objective coloring: ggnewscale", call. = FALSE)
  }
}

plot_frame <- function(df, coord_names) {
  out <- df
  out$.embedding_x <- suppressWarnings(as.numeric(out[[coord_names[[1L]]]]))
  out$.embedding_y <- suppressWarnings(as.numeric(out[[coord_names[[2L]]]]))
  out
}

square_embedding_limits <- function(plot_data, pad_frac = 0.05, coord_names = c("UMAP1", "UMAP2")) {
  coord_names <- as.character(coord_names)
  if (!all(coord_names %in% names(plot_data))) {
    stop("Plot data is missing coordinate columns: ", paste(setdiff(coord_names, names(plot_data)), collapse = ", "), call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(plot_data[[coord_names[[1L]]]]))
  y <- suppressWarnings(as.numeric(plot_data[[coord_names[[2L]]]]))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("Plot data must contain finite embedding coordinates.", call. = FALSE)
  x_range <- range(x)
  y_range <- range(y)
  span <- max(diff(x_range), diff(y_range))
  if (!is.finite(span) || span <= 0) span <- 1
  pad_frac <- suppressWarnings(as.numeric(pad_frac))
  if (!is.finite(pad_frac) || pad_frac < 0) pad_frac <- 0.05
  span <- span * (1 + 2 * pad_frac)
  x_center <- mean(x_range)
  y_center <- mean(y_range)
  list(
    xlim = c(x_center - span / 2, x_center + span / 2),
    ylim = c(y_center - span / 2, y_center + span / 2)
  )
}

add_plot_label_margin_rows <- function(plot_data, coord_names = c("UMAP1", "UMAP2"), pad_frac = 0.14) {
  coord_names <- as.character(coord_names)
  pad_frac <- suppressWarnings(as.numeric(pad_frac))
  if (!is.finite(pad_frac) || pad_frac < 0) pad_frac <- 0.14
  lims <- square_embedding_limits(plot_data, pad_frac = pad_frac, coord_names = coord_names)
  margin_rows <- plot_data[rep(1L, 4L), , drop = FALSE]
  for (nm in names(margin_rows)) margin_rows[[nm]] <- NA
  margin_rows[[coord_names[[1L]]]] <- c(lims$xlim[[1L]], lims$xlim[[2L]], lims$xlim[[1L]], lims$xlim[[2L]])
  margin_rows[[coord_names[[2L]]]] <- c(lims$ylim[[1L]], lims$ylim[[1L]], lims$ylim[[2L]], lims$ylim[[2L]])
  if ("point_type" %in% names(margin_rows)) margin_rows$point_type <- "plot_margin"
  if ("dataset" %in% names(margin_rows)) margin_rows$dataset <- "plot_margin"
  rbind(plot_data, margin_rows)
}

cluster_palette <- function(ids) {
  ids <- sort(unique(as.character(ids)))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (!length(ids)) return(stats::setNames(character(), character()))
  values <- grDevices::hcl.colors(length(ids), palette = "Dark 3")
  stats::setNames(values, ids)
}

cluster_hull_data <- function(df, cluster_col, expand = 1.035) {
  df <- df[!is.na(df[[cluster_col]]) & nzchar(as.character(df[[cluster_col]])), , drop = FALSE]
  if (!nrow(df)) return(data.frame())
  all_span <- max(
    diff(range(df$.embedding_x, finite = TRUE)),
    diff(range(df$.embedding_y, finite = TRUE))
  )
  point_radius <- if (is.finite(all_span) && all_span > 0) all_span * 0.012 else 0.1
  pieces <- lapply(split(df, as.character(df[[cluster_col]])), function(d) {
    xy <- unique(d[, c(".embedding_x", ".embedding_y"), drop = FALSE])
    if (nrow(xy) >= 3L) {
      idx <- grDevices::chull(xy$.embedding_x, xy$.embedding_y)
      h <- xy[c(idx, idx[[1L]]), , drop = FALSE]
      center <- colMeans(h[, c(".embedding_x", ".embedding_y"), drop = FALSE])
      h$.embedding_x <- center[[1L]] + (h$.embedding_x - center[[1L]]) * expand
      h$.embedding_y <- center[[2L]] + (h$.embedding_y - center[[2L]]) * expand
    } else if (nrow(xy) >= 1L) {
      theta <- seq(0, 2 * pi, length.out = 33L)
      h <- data.frame(
        .embedding_x = xy$.embedding_x[[1L]] + point_radius * cos(theta),
        .embedding_y = xy$.embedding_y[[1L]] + point_radius * sin(theta),
        stringsAsFactors = FALSE
      )
    } else {
      return(NULL)
    }
    h$cluster_id <- as.character(d[[cluster_col]][[1L]])
    h
  })
  rbind_fill_plain(Filter(Negate(is.null), pieces))
}

external_cluster_label_data <- function(clustered_plot_data,
                                        bounds_data,
                                        coord_names = c("UMAP1", "UMAP2"),
                                        cluster_col = "cluster_id",
                                        offset_frac = 0.07) {
  coord_names <- as.character(coord_names)
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  all_span <- max(diff(lims$xlim), diff(lims$ylim))
  if (!is.finite(all_span) || all_span <= 0) all_span <- 1
  margin <- all_span * 0.035
  offset <- all_span * offset_frac
  bx <- suppressWarnings(as.numeric(bounds_data[[coord_names[[1L]]]]))
  by <- suppressWarnings(as.numeric(bounds_data[[coord_names[[2L]]]]))
  global_center <- c(stats::median(bx, na.rm = TRUE), stats::median(by, na.rm = TRUE))
  if (any(!is.finite(global_center))) global_center <- c(mean(lims$xlim), mean(lims$ylim))
  cluster_ids <- sort(unique(as.character(clustered_plot_data[[cluster_col]])))
  cluster_ids <- cluster_ids[nzchar(cluster_ids) & !is.na(cluster_ids)]
  labels <- lapply(cluster_ids, function(cluster_id) {
    d <- clustered_plot_data[as.character(clustered_plot_data[[cluster_col]]) == cluster_id, coord_names, drop = FALSE]
    names(d) <- c(".embedding_x", ".embedding_y")
    d$.embedding_x <- suppressWarnings(as.numeric(d$.embedding_x))
    d$.embedding_y <- suppressWarnings(as.numeric(d$.embedding_y))
    d <- d[is.finite(d$.embedding_x) & is.finite(d$.embedding_y), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    center <- c(stats::median(d$.embedding_x), stats::median(d$.embedding_y))
    direction <- center - global_center
    norm <- sqrt(sum(direction^2))
    if (!is.finite(norm) || norm <= 1e-9) {
      far_i <- which.max((d$.embedding_x - center[[1L]])^2 + (d$.embedding_y - center[[2L]])^2)
      direction <- c(d$.embedding_x[[far_i]], d$.embedding_y[[far_i]]) - center
      norm <- sqrt(sum(direction^2))
    }
    if (!is.finite(norm) || norm <= 1e-9) {
      direction <- c(1, 0)
      norm <- 1
    }
    direction <- direction / norm
    boundary_i <- which.max((d$.embedding_x - center[[1L]]) * direction[[1L]] + (d$.embedding_y - center[[2L]]) * direction[[2L]])
    target_pos <- c(d$.embedding_x[[boundary_i]], d$.embedding_y[[boundary_i]])
    label_pos <- target_pos + direction * offset
    label_pos[[1L]] <- min(max(label_pos[[1L]], lims$xlim[[1L]] + margin), lims$xlim[[2L]] - margin)
    label_pos[[2L]] <- min(max(label_pos[[2L]], lims$ylim[[1L]] + margin), lims$ylim[[2L]] - margin)
    data.frame(
      cluster_id = cluster_id,
      .embedding_x = label_pos[[1L]],
      .embedding_y = label_pos[[2L]],
      .target_x = target_pos[[1L]],
      .target_y = target_pos[[2L]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, Filter(Negate(is.null), labels))
  if (is.null(out)) out <- data.frame()
  row.names(out) <- NULL
  out
}

add_cluster_outline_layers_external_labels <- function(plot,
                                                       clustered_plot_data,
                                                       bounds_data,
                                                       coord_names = c("UMAP1", "UMAP2"),
                                                       cluster_col = "cluster_id",
                                                       linewidth = 0.72) {
  clustered_plot_data <- plot_frame(clustered_plot_data, coord_names = coord_names)
  hulls <- cluster_hull_data(clustered_plot_data, cluster_col = cluster_col)
  if (!nrow(hulls)) return(plot)
  pal <- cluster_palette(hulls$cluster_id)
  for (cluster_id in names(pal)) {
    h <- hulls[hulls$cluster_id == cluster_id, , drop = FALSE]
    plot <- plot +
      ggplot2::geom_path(
        data = h,
        ggplot2::aes(x = .embedding_x, y = .embedding_y),
        inherit.aes = FALSE,
        color = pal[[cluster_id]],
        linewidth = linewidth,
        linetype = "dashed",
        lineend = "round",
        show.legend = FALSE
      )
  }
  labels <- external_cluster_label_data(
    clustered_plot_data = clustered_plot_data,
    bounds_data = bounds_data,
    coord_names = coord_names,
    cluster_col = cluster_col
  )
  if (!nrow(labels)) return(plot)
  labels$.nudge_x <- labels$.embedding_x - labels$.target_x
  labels$.nudge_y <- labels$.embedding_y - labels$.target_y
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    plot <- plot + ggnewscale::new_scale_color()
  }
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(
      plot +
        ggrepel::geom_label_repel(
          data = labels,
          ggplot2::aes(x = .target_x, y = .target_y, label = cluster_id, color = cluster_id),
          inherit.aes = FALSE,
          fill = "white",
          label.size = 0.18,
          size = 2.4,
          fontface = "bold",
          box.padding = 0.4,
          point.padding = 0.15,
          min.segment.length = 0,
          segment.size = 0.25,
          nudge_x = labels$.nudge_x,
          nudge_y = labels$.nudge_y,
          seed = 123,
          max.overlaps = Inf,
          xlim = lims$xlim,
          ylim = lims$ylim,
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = pal, guide = "none")
    )
  }
  plot +
    ggplot2::geom_segment(
      data = labels,
      ggplot2::aes(x = .target_x, y = .target_y, xend = .embedding_x, yend = .embedding_y, color = cluster_id),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.85,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = labels,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, label = cluster_id, color = cluster_id),
      inherit.aes = FALSE,
      fill = "white",
      linewidth = 0.18,
      size = 2.4,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = pal, guide = "none")
}

base_embedding_plot <- function(plot_data,
                                coord_names,
                                axis_labels,
                                initial_size = 0.22,
                                best_size = 1.25) {
  require_pooled_objective_plotting()
  lims <- square_embedding_limits(plot_data, coord_names = coord_names)
  dat <- plot_frame(plot_data, coord_names = coord_names)
  dat$dataset <- factor(as.character(dat$dataset), levels = c("invivo", "invitro", "plot_margin"))
  dat$point_type <- as.character(dat$point_type)
  initial <- dat[dat$point_type == "initial", , drop = FALSE]
  initial_invivo <- initial[initial$dataset == "invivo", , drop = FALSE]
  initial_invitro <- initial[initial$dataset == "invitro", , drop = FALSE]
  best_invivo <- dat[dat$point_type == "best" & dat$dataset == "invivo", , drop = FALSE]
  best_invitro <- dat[dat$point_type == "best" & dat$dataset == "invitro", , drop = FALSE]
  initial_invivo_color <- "#7FA2FF"
  initial_invitro_color <- "#FF4FBF"
  initial_invivo_alpha <- 0.48
  initial_invitro_alpha <- 0.45
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = initial_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 16,
      color = initial_invivo_color,
      alpha = initial_invivo_alpha,
      size = initial_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = initial_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y),
      shape = 17,
      color = initial_invitro_color,
      alpha = initial_invitro_alpha,
      size = initial_size,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = "Initial samples",
      values = c(invivo = 16, invitro = 17),
      labels = c(invivo = "in vivo", invitro = "in vitro")
    ) +
    ggplot2::geom_point(
      data = best_invivo,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective, shape = dataset),
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vivo\nobjective",
      low = "#2C7BB6",
      high = "#FDE725",
      guide = ggplot2::guide_colorbar(order = 2, barheight = ggplot2::unit(26, "mm"))
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = best_invitro,
      ggplot2::aes(x = .embedding_x, y = .embedding_y, color = objective, shape = dataset),
      alpha = 0.95,
      size = best_size,
      stroke = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_gradient(
      name = "in vitro\nobjective",
      low = "#1A9850",
      high = "#D73027",
      guide = ggplot2::guide_colorbar(order = 3, barheight = ggplot2::unit(26, "mm"))
    )
  p +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]]) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          color = c(initial_invivo_color, initial_invitro_color),
          alpha = c(initial_invivo_alpha, initial_invitro_alpha),
          size = 3
        )
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      legend.box = "vertical",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

subcluster_embedding_plot <- function(subclustered_best, coord_names, axis_labels, title = NULL) {
  require_plotting()
  dat <- plot_frame(subclustered_best, coord_names = coord_names)
  bounds_data <- add_plot_label_margin_rows(dat, coord_names = coord_names, pad_frac = 0.2)
  lims <- square_embedding_limits(bounds_data, coord_names = coord_names)
  dat$subcluster_id <- as.character(dat$subcluster_id)
  dat$dataset <- factor(as.character(dat$dataset), levels = c("invivo", "invitro"))
  pal <- cluster_palette(dat$subcluster_id)
  ggplot2::ggplot(dat, ggplot2::aes(x = .embedding_x, y = .embedding_y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = subcluster_id, shape = dataset),
      size = 1.7,
      alpha = 0.95,
      stroke = 0.15
    ) +
    ggplot2::scale_color_manual(name = "subcluster", values = pal, guide = "none") +
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c(invivo = 16, invitro = 17),
      labels = c(invivo = "in vivo", invitro = "in vitro")
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = axis_labels[[1L]], y = axis_labels[[2L]], title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

save_plot_pair <- function(plot, stem, width = 7.4, height = 6.5) {
  require_plotting()
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 220, bg = "white")
  invisible(stem)
}

write_cluster_figures <- function(full_marked,
                                  clustered_best,
                                  subclustered_best,
                                  figure_dir,
                                  output_stem,
                                  reduction,
                                  coord_names) {
  require_plotting()
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  axis_labels <- reduction_axis_labels(reduction)
  plot_bounds <- add_plot_label_margin_rows(full_marked, coord_names = coord_names, pad_frac = 0.14)
  main_plot <- base_embedding_plot(
    plot_bounds,
    coord_names = coord_names,
    axis_labels = axis_labels
  )
  main_plot <- add_cluster_outline_layers_external_labels(
    main_plot,
    clustered_best,
    bounds_data = plot_bounds,
    coord_names = coord_names,
    cluster_col = "cluster_id"
  )
  save_plot_pair(main_plot, file.path(figure_dir, output_stem), width = 7.4, height = 6.5)

  sub_dir <- file.path(figure_dir, "Subclusters")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  split_vals <- sort(unique(as.character(subclustered_best$primary_cluster_id)))
  for (primary_id in split_vals) {
    d <- subclustered_best[as.character(subclustered_best$primary_cluster_id) == primary_id, , drop = FALSE]
    if (!nrow(d)) next
    sub_bounds <- add_plot_label_margin_rows(d, coord_names = coord_names, pad_frac = 0.2)
    p <- subcluster_embedding_plot(
      d,
      coord_names = coord_names,
      axis_labels = axis_labels,
      title = paste0(unique(d$method)[[1L]], " ", primary_id, " subclusters")
    )
    p <- add_cluster_outline_layers_external_labels(
      p,
      d,
      bounds_data = sub_bounds,
      coord_names = coord_names,
      cluster_col = "subcluster_id",
      linewidth = 0.72
    )
    save_plot_pair(
      p,
      file.path(sub_dir, paste0(reduction_file_suffix(reduction), "_", primary_id, "_subclusters")),
      width = 5.4,
      height = 4.8
    )
  }
  invisible(TRUE)
}

analyze_embedding <- function(reduction,
                              coordinate_csv,
                              feature_data,
                              output_dir,
                              cluster_seed,
                              cluster_k_min,
                              cluster_k_max,
                              silhouette_sample_n,
                              subcluster_seed,
                              subcluster_k_min,
                              subcluster_k_max,
                              subcluster_min_n) {
  reduction <- normalize_reduction(reduction)
  coord_names <- reduction_coordinate_names(reduction)
  plot_data <- read_csv_plain(coordinate_csv)
  required <- c(coord_names, "dataset", "point_type", "seed", "objective")
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) stop("Coordinate CSV is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  plot_data$dataset <- as.character(plot_data$dataset)
  plot_data$point_type <- as.character(plot_data$point_type)
  plot_data$seed <- suppressWarnings(as.integer(plot_data$seed))
  plot_data$objective <- suppressWarnings(as.numeric(plot_data$objective))
  for (nm in coord_names) plot_data[[nm]] <- suppressWarnings(as.numeric(plot_data[[nm]]))
  if (!all(stats::complete.cases(plot_data[, coord_names, drop = FALSE]))) {
    stop("Coordinate CSV contains non-finite ", reduction, " coordinates.", call. = FALSE)
  }
  specs <- cluster_dataset_specs()
  cluster_results <- lapply(seq_len(nrow(specs)), function(i) {
    cluster_best_dataset(
      plot_data = plot_data,
      dataset = specs$dataset[[i]],
      dataset_label = specs$dataset_label[[i]],
      cluster_prefix = specs$cluster_prefix[[i]],
      coord_names = coord_names,
      cluster_seed = as.integer(cluster_seed) + i - 1L,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      silhouette_sample_n = silhouette_sample_n
    )
  })
  names(cluster_results) <- specs$dataset
  clustered_best <- do.call(rbind, lapply(cluster_results, `[[`, "best"))
  silhouette_all <- do.call(rbind, lapply(cluster_results, `[[`, "silhouette"))
  cluster_summary <- do.call(rbind, lapply(cluster_results, `[[`, "summary"))

  full_marked <- plot_data
  full_marked$dataset_label <- NA_character_
  full_marked$cluster_scope <- NA_character_
  full_marked$cluster_source <- NA_character_
  full_marked$cluster_prefix <- NA_character_
  full_marked$cluster_base_id <- NA_character_
  full_marked$cluster_id <- NA_character_
  full_marked$cluster_num <- NA_integer_
  full_marked$cluster_k <- NA_integer_
  full_marked$cluster_silhouette_avg <- NA_real_
  full_marked$cluster_silhouette_sample_n <- NA_integer_
  cluster_cols <- c(
    "dataset_label", "cluster_scope", "cluster_source", "cluster_prefix",
    "cluster_base_id", "cluster_id", "cluster_num", "cluster_k",
    "cluster_silhouette_avg", "cluster_silhouette_sample_n"
  )
  for (cluster_result in cluster_results) {
    full_marked[cluster_result$best_idx, cluster_cols] <- cluster_result$best[, cluster_cols, drop = FALSE]
  }

  seed_groups <- best_seed_group_table(clustered_best, reduction = reduction, coord_names = coord_names)
  subclusters <- analyze_best_subclusters(
    clustered_best = clustered_best,
    feature_data = feature_data,
    reduction = reduction,
    coord_names = coord_names,
    subcluster_seed = subcluster_seed,
    subcluster_k_min = subcluster_k_min,
    subcluster_k_max = subcluster_k_max,
    subcluster_min_n = subcluster_min_n,
    silhouette_sample_n = silhouette_sample_n
  )
  table_dir <- file.path(output_dir, "Tables")
  figure_dir <- file.path(output_dir, "Figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- reduction_file_suffix(reduction)
  stem <- paste0("pooled_invivo_invitro_initial_vs_best_", suffix, "_best_clusters")
  write_csv(full_marked, file.path(table_dir, paste0(stem, "_full_coordinates.csv")))
  write_csv(clustered_best, file.path(table_dir, paste0(stem, "_best_coordinates.csv")))
  write_csv(silhouette_all, file.path(table_dir, paste0(stem, "_silhouette.csv")))
  write_csv(cluster_summary, file.path(table_dir, paste0(stem, "_cluster_summary.csv")))
  write_csv(seed_groups, file.path(table_dir, paste0(stem, "_best_seed_groups.csv")))
  write_csv(
    data.frame(
      key = c(
        "reduction", "coordinate_csv", "output_dir", "n_full_rows",
        "n_invivo_best_rows", "n_invitro_best_rows", "coordinate_columns",
        "cluster_seed", "cluster_k_min", "cluster_k_max", "invivo_selected_k",
        "invitro_selected_k", "invivo_selected_average_silhouette",
        "invitro_selected_average_silhouette", "subcluster_seed",
        "subcluster_k_min", "subcluster_k_max", "subcluster_min_n", "created_at"
      ),
      value = c(
        reduction, coordinate_csv, output_dir, nrow(plot_data),
        nrow(cluster_results$invivo$best), nrow(cluster_results$invitro$best),
        paste(coord_names, collapse = ","),
        cluster_seed, cluster_k_min, cluster_k_max,
        cluster_results$invivo$selected_summary$k[[1L]],
        cluster_results$invitro$selected_summary$k[[1L]],
        cluster_results$invivo$selected_summary$average_silhouette[[1L]],
        cluster_results$invitro$selected_summary$average_silhouette[[1L]],
        subcluster_seed, subcluster_k_min, subcluster_k_max, subcluster_min_n,
        format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      ),
      stringsAsFactors = FALSE
    ),
    file.path(table_dir, paste0(stem, "_metadata.csv"))
  )
  write_csv(subclusters$best, file.path(table_dir, paste0(stem, "_best_subclusters.csv")))
  write_csv(subclusters$seed_groups, file.path(table_dir, paste0(stem, "_best_subcluster_seed_groups.csv")))
  write_csv(subclusters$summary, file.path(table_dir, paste0(stem, "_best_subcluster_summary.csv")))
  write_csv(subclusters$silhouette, file.path(table_dir, paste0(stem, "_best_subcluster_silhouette.csv")))
  write_csv(subclusters$feature_metadata, file.path(table_dir, paste0(stem, "_best_subcluster_feature_metadata.csv")))
  write_csv(subclusters$zscore_features, file.path(table_dir, paste0(stem, "_best_subcluster_zscore_features.csv")))
  write_cluster_figures(
    full_marked = full_marked,
    clustered_best = clustered_best,
    subclustered_best = subclusters$best,
    figure_dir = figure_dir,
    output_stem = stem,
    reduction = reduction,
    coord_names = coord_names
  )
  list(seed_groups = seed_groups, subcluster_seed_groups = subclusters$seed_groups, subcluster_summary = subclusters$summary)
}

seed_dir_lookup <- function(seed_dirs, dataset) {
  data.frame(
    dataset = dataset,
    seed = vapply(seed_dirs, seed_from_dir, integer(1L)),
    seed_dir = normalizePath(seed_dirs, mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

selected_representatives <- function(summary_df, seed_dirs_df) {
  rows <- summary_df[order(summary_df$method, summary_df$dataset, summary_df$primary_cluster_num, summary_df$subcluster_num), , drop = FALSE]
  rows$seed <- as.integer(rows$objective_min_seed)
  key <- seed_key(rows$dataset, rows$seed)
  seed_dirs_df$.key <- seed_key(seed_dirs_df$dataset, seed_dirs_df$seed)
  idx <- match(key, seed_dirs_df$.key)
  if (any(is.na(idx))) stop("Could not map representative seeds to seed directories: ", paste(head(key[is.na(idx)], 20), collapse = ", "), call. = FALSE)
  rows$seed_dir <- seed_dirs_df$seed_dir[idx]
  rows$representative_rank <- ave(seq_len(nrow(rows)), rows$method, rows$dataset, FUN = seq_along)
  rows
}

global_invitro_best_anchor <- function(invitro_best_csv, seed_dirs_df) {
  best <- read_csv_plain(invitro_best_csv)
  if (!all(c("seed", "objective") %in% names(best))) {
    stop("In vitro best-parameter table must contain seed and objective columns: ", invitro_best_csv, call. = FALSE)
  }
  best$seed <- suppressWarnings(as.integer(best$seed))
  best$objective <- suppressWarnings(as.numeric(best$objective))
  finite <- is.finite(best$objective) & !is.na(best$seed)
  if (!any(finite)) stop("No finite in vitro objective values found in: ", invitro_best_csv, call. = FALSE)
  best <- best[finite, , drop = FALSE]
  best <- best[order(best$objective, best$seed), , drop = FALSE]
  seed <- as.integer(best$seed[[1L]])
  seed_dirs_df$.key <- seed_key(seed_dirs_df$dataset, seed_dirs_df$seed)
  idx <- match(seed_key("invitro", seed), seed_dirs_df$.key)
  if (is.na(idx)) stop("Could not map global best in vitro seed to a seed directory: ", seed, call. = FALSE)
  data.frame(
    dataset = "invitro",
    seed = seed,
    seed_dir = as.character(seed_dirs_df$seed_dir[[idx]]),
    objective = as.numeric(best$objective[[1L]]),
    representative_rank = 1L,
    family = "global_best",
    stringsAsFactors = FALSE
  )
}

manifest_columns <- c(
  "warmup_label", "phase", "invivo_family", "invivo_rank", "invivo_seed", "invivo_seed_dir",
  "invitro_family", "invitro_rank", "invitro_seed", "invitro_seed_dir",
  "selection_reason", "joint_run_prefix", "joint_soft_coupling_parameters_table"
)

normalize_pairing_policy <- function(pairing_policy) {
  pairing_policy <- tolower(trimws(as.character(pairing_policy)))
  pairing_policy <- gsub("-", "_", pairing_policy, fixed = TRUE)
  aliases <- c(
    invitro_best_to_invivo_subcluster = "invitro_best_to_invivo_subclusters",
    invivo_subclusters_to_invitro_best = "invitro_best_to_invivo_subclusters",
    invitro_best_by_invivo_subclusters = "invitro_best_to_invivo_subclusters"
  )
  if (pairing_policy %in% names(aliases)) pairing_policy <- aliases[[pairing_policy]]
  pairing_policy
}

row_scalar <- function(row, name, default = "") {
  if (!name %in% names(row)) return(default)
  val <- row[[name]][[1L]]
  if (is.null(val) || length(val) == 0L || is.na(val)) default else as.character(val)
}

cluster_subcluster_token <- function(row) {
  primary <- row_scalar(row, "primary_cluster_base_id")
  sub <- row_scalar(row, "subcluster_base_id")
  token <- paste0(primary, sub)
  if (!nzchar(token)) token <- row_scalar(row, "subcluster_id", "cluster")
  sanitize_label(token, fallback = "cluster")
}

warmup_label_for_pair <- function(method, invivo_row, invitro_row, include_invitro_cluster = FALSE) {
  invivo_seed <- as.integer(invivo_row$seed[[1L]])
  invitro_seed <- as.integer(invitro_row$seed[[1L]])
  invivo_token <- paste0("vi_seed", invivo_seed, "_", cluster_subcluster_token(invivo_row))
  invitro_token <- paste0("vt_seed", invitro_seed)
  if (isTRUE(include_invitro_cluster)) {
    invitro_token <- paste0(invitro_token, "_", cluster_subcluster_token(invitro_row))
  }
  sanitize_label(paste(method, invivo_token, invitro_token, sep = "_"))
}

empty_manifest <- function() {
  as.data.frame(setNames(rep(list(character()), length(manifest_columns)), manifest_columns), stringsAsFactors = FALSE)[0L, , drop = FALSE]
}

manifest_row <- function(warmup_label,
                         method,
                         invivo_row,
                         invitro_row,
                         invitro_family,
                         invitro_rank,
                         selection_reason,
                         out_dir) {
  data.frame(
    warmup_label = warmup_label,
    phase = method,
    invivo_family = as.character(invivo_row$subcluster_id[[1L]]),
    invivo_rank = as.integer(invivo_row$representative_rank[[1L]]),
    invivo_seed = as.integer(invivo_row$seed[[1L]]),
    invivo_seed_dir = as.character(invivo_row$seed_dir[[1L]]),
    invitro_family = invitro_family,
    invitro_rank = as.integer(invitro_rank),
    invitro_seed = as.integer(invitro_row$seed[[1L]]),
    invitro_seed_dir = as.character(invitro_row$seed_dir[[1L]]),
    selection_reason = selection_reason,
    joint_run_prefix = paste0("fit_joint_", warmup_label),
    joint_soft_coupling_parameters_table = file.path(out_dir, "joint_soft_coupling_tables", paste0("joint_soft_coupling_parameters_table__", warmup_label, ".csv")),
    stringsAsFactors = FALSE
  )
}

build_manifest_from_representatives <- function(reps,
                                                out_dir,
                                                pairing_policy = "cartesian_by_method",
                                                deduplicate_pairs = FALSE,
                                                invitro_best_anchor = NULL) {
  pairing_policy <- normalize_pairing_policy(pairing_policy)
  if (!pairing_policy %in% c("cartesian_by_method", "invitro_best_to_invivo_subclusters")) {
    stop(
      "--pairing_policy must be cartesian_by_method or invitro_best_to_invivo_subclusters, got: ",
      pairing_policy,
      call. = FALSE
    )
  }
  methods <- unique(as.character(reps$method))
  rows <- list()
  idx <- 1L
  for (method in methods) {
    invivo <- reps[reps$method == method & reps$dataset == "invivo", , drop = FALSE]
    if (!nrow(invivo)) next
    if (identical(pairing_policy, "invitro_best_to_invivo_subclusters")) {
      if (is.null(invitro_best_anchor) || !nrow(invitro_best_anchor)) {
        stop("pairing_policy=invitro_best_to_invivo_subclusters requires a global in vitro best anchor.", call. = FALSE)
      }
      invitro <- invitro_best_anchor[1L, , drop = FALSE]
      for (i in seq_len(nrow(invivo))) {
        warmup_label <- warmup_label_for_pair(
          method = method,
          invivo_row = invivo[i, , drop = FALSE],
          invitro_row = invitro,
          include_invitro_cluster = FALSE
        )
        rows[[idx]] <- manifest_row(
          warmup_label = warmup_label,
          method = method,
          invivo_row = invivo[i, , drop = FALSE],
          invitro_row = invitro,
          invitro_family = "global_best",
          invitro_rank = 1L,
          selection_reason = "landscape_invivo_subcluster_objective_min_seed_to_global_invitro_objective_min_seed",
          out_dir = out_dir
        )
        idx <- idx + 1L
      }
      next
    }

    invitro <- reps[reps$method == method & reps$dataset == "invitro", , drop = FALSE]
    if (!nrow(invitro)) next
    for (i in seq_len(nrow(invivo))) {
      for (j in seq_len(nrow(invitro))) {
        warmup_label <- warmup_label_for_pair(
          method = method,
          invivo_row = invivo[i, , drop = FALSE],
          invitro_row = invitro[j, , drop = FALSE],
          include_invitro_cluster = TRUE
        )
        rows[[idx]] <- manifest_row(
          warmup_label = warmup_label,
          method = method,
          invivo_row = invivo[i, , drop = FALSE],
          invitro_row = invitro[j, , drop = FALSE],
          invitro_family = as.character(invitro$subcluster_id[[j]]),
          invitro_rank = as.integer(invitro$representative_rank[[j]]),
          selection_reason = "landscape_second_level_cluster_objective_min_seed",
          out_dir = out_dir
        )
        idx <- idx + 1L
      }
    }
  }
  if (!length(rows)) {
    return(empty_manifest())
  }
  manifest <- do.call(rbind, rows)
  if (isTRUE(deduplicate_pairs)) {
    key <- paste(manifest$invivo_seed, manifest$invitro_seed, sep = "::")
    manifest <- manifest[!duplicated(key), , drop = FALSE]
    row.names(manifest) <- NULL
  }
  manifest[, manifest_columns, drop = FALSE]
}

landscape_analysis_root <- function(out_dir) {
  file.path(out_dir, "landscape_subcluster")
}

landscape_cluster_output_dir <- function(analysis_root) {
  file.path(analysis_root, "pooled_invivo_invitro", "full_data_in_vivo_clustring")
}

landscape_combined_table_dir <- function(analysis_root) {
  file.path(landscape_cluster_output_dir(analysis_root), "Tables")
}

warmup_cross_validation_root <- function(out_dir) {
  file.path(out_dir, "cross_validation", "best_fit_parameter_feature")
}

warmup_dense_grid_dir <- function(out_dir) {
  file.path(
    warmup_cross_validation_root(out_dir),
    "03_dense-grid_monotonicity_classification", "monotonicity_classification",
    "dense-grid_monotonicity_classification"
  )
}

landscape_input_tables_path <- function(out_dir) {
  file.path(out_dir, "landscape_subcluster_input_tables.tsv")
}

write_landscape_input_tables <- function(out_dir,
                                         invivo_tables,
                                         invitro_tables,
                                         invivo_run_dir,
                                         invitro_run_dir,
                                         analysis_root) {
  write_tsv(
    data.frame(
      key = c(
        "analysis_root", "seed_table_dir",
        "invivo_run_dir", "invitro_run_dir",
        "invivo_best_csv", "invivo_initial_csv",
        "invitro_best_csv", "invitro_initial_csv"
      ),
      value = c(
        analysis_root, file.path(analysis_root, "SeedParameterTables"),
        invivo_run_dir, invitro_run_dir,
        invivo_tables$best_csv, invivo_tables$initial_csv,
        invitro_tables$best_csv, invitro_tables$initial_csv
      ),
      stringsAsFactors = FALSE
    ),
    landscape_input_tables_path(out_dir)
  )
}

read_key_value_table <- function(path) {
  if (!file.exists(path)) return(stats::setNames(character(), character()))
  tab <- read_tsv(path)
  if (!all(c("key", "value") %in% names(tab))) return(stats::setNames(character(), character()))
  stats::setNames(as.character(tab$value), as.character(tab$key))
}

table_value <- function(tab, key, default = "") {
  if (key %in% names(tab) && nzchar(tab[[key]])) tab[[key]] else default
}

curve_class_by_seed_path <- function(out_dir) {
  file.path(warmup_dense_grid_dir(out_dir), "tables", "fixed_o2_ploidy_monotonicity_by_seed.tsv")
}

read_curve_class_table <- function(path, class_col = "curve_class") {
  tab <- read_tsv(path)
  if (!"seed_number" %in% names(tab)) {
    if (!"seed_id" %in% names(tab)) stop("Curve-class table must contain seed_number or seed_id: ", path, call. = FALSE)
    tab$seed_number <- suppressWarnings(as.integer(sub("^seed", "", as.character(tab$seed_id))))
  }
  if (!class_col %in% names(tab)) stop("Curve-class table is missing column ", class_col, ": ", path, call. = FALSE)
  tab$seed_number <- suppressWarnings(as.integer(tab$seed_number))
  tab[[class_col]] <- trimws(as.character(tab[[class_col]]))
  tab <- tab[is.finite(tab$seed_number) & nzchar(tab[[class_col]]), , drop = FALSE]
  tab <- tab[!duplicated(tab$seed_number), , drop = FALSE]
  tab[, unique(c("seed_number", class_col, intersect(c("seed_id", "final_interpretation_class", "monotonicity_reliability"), names(tab)))), drop = FALSE]
}

filter_invivo_summary_by_curve_class <- function(summary_df,
                                                 seed_groups_df,
                                                 curve_class_table,
                                                 target_class = "monotone_increasing",
                                                 class_col = "curve_class") {
  curve_class_table$seed_number <- suppressWarnings(as.integer(curve_class_table$seed_number))
  curve_class_table[[class_col]] <- as.character(curve_class_table[[class_col]])
  rows <- list()
  audit <- list()
  out_i <- 0L
  audit_i <- 0L

  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, , drop = FALSE]
    if (!identical(as.character(row$dataset[[1L]]), "invivo")) {
      row$representative_selection_rule <- "objective_min_seed"
      row$representative_curve_class <- NA_character_
      row$n_curve_class_candidate_seeds <- NA_integer_
      row$original_objective_min_seed <- row$objective_min_seed
      row$original_objective_min <- row$objective_min
      out_i <- out_i + 1L
      rows[[out_i]] <- row
      next
    }

    d <- seed_groups_df[
      as.character(seed_groups_df$method) == as.character(row$method[[1L]]) &
        as.character(seed_groups_df$dataset) == "invivo" &
        as.character(seed_groups_df$subcluster_id) == as.character(row$subcluster_id[[1L]]),
      ,
      drop = FALSE
    ]
    if (!nrow(d)) {
      stop("Could not find seed-group rows for in vivo subcluster: ", row$method[[1L]], " ", row$subcluster_id[[1L]], call. = FALSE)
    }
    d$seed <- suppressWarnings(as.integer(d$seed))
    d$objective <- suppressWarnings(as.numeric(d$objective))
    idx <- match(d$seed, curve_class_table$seed_number)
    d$curve_class_for_filter <- curve_class_table[[class_col]][idx]
    keep <- d[!is.na(d$curve_class_for_filter) & d$curve_class_for_filter == target_class, , drop = FALSE]
    audit_i <- audit_i + 1L
    audit[[audit_i]] <- data.frame(
      method = as.character(row$method[[1L]]),
      dataset = "invivo",
      primary_cluster_id = as.character(row$primary_cluster_id[[1L]]),
      subcluster_id = as.character(row$subcluster_id[[1L]]),
      n_seeds = nrow(d),
      target_curve_class = target_class,
      n_curve_class_candidate_seeds = nrow(keep),
      original_objective_min_seed = as.integer(row$objective_min_seed[[1L]]),
      original_objective_min = suppressWarnings(as.numeric(row$objective_min[[1L]])),
      selected_seed = if (nrow(keep)) keep$seed[order(keep$objective, keep$seed)][[1L]] else NA_integer_,
      selected_objective = if (nrow(keep)) keep$objective[order(keep$objective, keep$seed)][[1L]] else NA_real_,
      status = if (nrow(keep)) "selected" else "skipped_no_matching_curve_class",
      stringsAsFactors = FALSE
    )
    if (!nrow(keep)) next
    keep <- keep[order(keep$objective, keep$seed), , drop = FALSE]
    selected <- keep[1L, , drop = FALSE]
    row$original_objective_min_seed <- row$objective_min_seed
    row$original_objective_min <- row$objective_min
    row$objective_min_seed <- as.integer(selected$seed[[1L]])
    row$objective_min <- suppressWarnings(as.numeric(selected$objective[[1L]]))
    row$representative_selection_rule <- paste0("objective_min_seed_with_curve_class_", target_class)
    row$representative_curve_class <- target_class
    row$n_curve_class_candidate_seeds <- nrow(keep)
    out_i <- out_i + 1L
    rows[[out_i]] <- row
  }

  filtered <- if (length(rows)) do.call(rbind, rows) else summary_df[0L, , drop = FALSE]
  audit_df <- if (length(audit)) do.call(rbind, audit) else data.frame()
  list(summary = filtered, audit = audit_df)
}

finalize_landscape_pairs_core <- function(out_dir,
                                          project_root,
                                          invivo_run_dir,
                                          invitro_run_dir,
                                          analysis_root,
                                          result_root,
                                          reductions,
                                          umap_seed,
                                          tsne_seed,
                                          cluster_seed,
                                          subcluster_seed,
                                          cluster_k_min,
                                          cluster_k_max,
                                          subcluster_k_min,
                                          subcluster_k_max,
                                          subcluster_min_n,
                                          pairing_policy,
                                          deduplicate_pairs,
                                          reference_dir,
                                          invivo_curve_filter = FALSE,
                                          invivo_curve_class = "monotone_increasing",
                                          invivo_curve_class_table = "",
                                          invitro_best_csv = "") {
  combined_table_dir <- landscape_combined_table_dir(analysis_root)
  summary_path <- file.path(combined_table_dir, "pooled_invivo_invitro_best_subcluster_summary_by_method.csv")
  seed_groups_path <- file.path(combined_table_dir, "pooled_invivo_invitro_best_subclusters_by_method.csv")
  if (!file.exists(summary_path)) stop("Missing subcluster summary table: ", summary_path, call. = FALSE)
  if (!file.exists(seed_groups_path)) stop("Missing subcluster seed-group table: ", seed_groups_path, call. = FALSE)
  all_subcluster_summary <- read_csv_plain(summary_path)
  all_subcluster_seed_groups <- read_csv_plain(seed_groups_path)

  seed_dirs_df <- rbind(
    seed_dir_lookup(list_seed_dirs(invivo_run_dir), "invivo"),
    seed_dir_lookup(list_seed_dirs(invitro_run_dir), "invitro")
  )
  selection_summary <- all_subcluster_summary
  if (isTRUE(invivo_curve_filter)) {
    if (!nzchar(invivo_curve_class_table)) {
      invivo_curve_class_table <- curve_class_by_seed_path(out_dir)
    }
    if (!file.exists(invivo_curve_class_table)) {
      stop("Missing in vivo curve-class table for representative filtering: ", invivo_curve_class_table, call. = FALSE)
    }
    curve_classes <- read_curve_class_table(invivo_curve_class_table, class_col = "curve_class")
    filtered <- filter_invivo_summary_by_curve_class(
      summary_df = selection_summary,
      seed_groups_df = all_subcluster_seed_groups,
      curve_class_table = curve_classes,
      target_class = invivo_curve_class,
      class_col = "curve_class"
    )
    selection_summary <- filtered$summary
    write_tsv(filtered$audit, file.path(out_dir, "landscape_subcluster_invivo_curve_class_representative_audit.tsv"))
    if (!any(as.character(selection_summary$dataset) == "invivo")) {
      stop("No in vivo subcluster representatives remained after curve-class filter: ", invivo_curve_class, call. = FALSE)
    }
  }

  reps <- selected_representatives(selection_summary, seed_dirs_df)
  write_tsv(reps, file.path(out_dir, "landscape_subcluster_selected_representatives.tsv"))
  if (!nzchar(invitro_best_csv)) {
    input_tables <- read_key_value_table(landscape_input_tables_path(out_dir))
    invitro_best_csv <- table_value(
      input_tables,
      "invitro_best_csv",
      file.path(analysis_root, "SeedParameterTables", "invitro_best_params_by_seed.csv")
    )
  }
  invitro_best_anchor <- global_invitro_best_anchor(invitro_best_csv, seed_dirs_df)
  write_tsv(invitro_best_anchor, file.path(out_dir, "landscape_subcluster_invitro_best_anchor.tsv"))
  manifest <- build_manifest_from_representatives(
    reps = reps,
    out_dir = out_dir,
    pairing_policy = pairing_policy,
    deduplicate_pairs = deduplicate_pairs,
    invitro_best_anchor = invitro_best_anchor
  )
  write_tsv(manifest, file.path(out_dir, "multi_warmup_manifest.tsv"))
  write_seed_plan_mode(
    out_dir,
    mode = if (identical(pairing_policy, "invitro_best_to_invivo_subclusters")) {
      "landscape_subcluster_invivo_to_invitro_best"
    } else {
      "landscape_subcluster_paired"
    },
    warmup_pairs = nrow(manifest),
    pairing_policy = pairing_policy,
    deduplicate_pairs = deduplicate_pairs
  )
  write_validation_summary(reference_dir, all_subcluster_summary, out_dir)
  write_tsv(
    data.frame(
      key = c(
        "pair_method", "invivo_run_dir", "invitro_run_dir", "analysis_root", "result_root",
        "reductions", "umap_seed", "tsne_seed", "cluster_seed", "subcluster_seed",
        "cluster_k_min", "cluster_k_max", "subcluster_k_min", "subcluster_k_max",
        "subcluster_min_n", "pairing_policy", "deduplicate_pairs",
        "invivo_curve_filter", "invivo_curve_class", "invivo_curve_class_table",
        "invitro_best_seed", "invitro_best_objective", "warmup_pairs"
      ),
      value = as.character(c(
        "landscape_subcluster", invivo_run_dir, invitro_run_dir, analysis_root, result_root,
        paste(reductions, collapse = ","), umap_seed, tsne_seed, cluster_seed, subcluster_seed,
        cluster_k_min, cluster_k_max, subcluster_k_min, subcluster_k_max,
        subcluster_min_n, pairing_policy, deduplicate_pairs,
        invivo_curve_filter, invivo_curve_class, invivo_curve_class_table,
        invitro_best_anchor$seed[[1L]], invitro_best_anchor$objective[[1L]], nrow(manifest)
      )),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "landscape_subcluster_pair_settings.tsv")
  )
  message("Wrote multi-warmup manifest: ", file.path(out_dir, "multi_warmup_manifest.tsv"))
  message("Selected representatives: ", nrow(reps), "; warm-up pairs: ", nrow(manifest))
  invisible(manifest)
}

write_seed_plan_mode <- function(out_dir, mode, warmup_pairs, pairing_policy, deduplicate_pairs) {
  write_tsv(
    data.frame(
      field = c("mode", "warmup_pairs", "pairing_policy", "deduplicate_pairs"),
      value = as.character(c(mode, warmup_pairs, pairing_policy, deduplicate_pairs)),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "multi_warmup_seed_plan_mode.tsv")
  )
}

write_validation_summary <- function(reference_dir, generated_summary, out_dir) {
  if (!nzchar(reference_dir)) return(invisible(NULL))
  reference_path <- file.path(reference_dir, "Tables", "pooled_invivo_invitro_best_subcluster_summary_by_method.csv")
  if (!file.exists(reference_path)) {
    warning("Reference subcluster summary not found: ", reference_path)
    return(invisible(NULL))
  }
  ref <- read_csv_plain(reference_path)
  key_cols <- c("method", "dataset", "primary_cluster_id", "subcluster_id")
  compare_cols <- c(key_cols, "objective_min_seed", "n_seeds", "subcluster_k")
  if (!all(compare_cols %in% names(ref)) || !all(compare_cols %in% names(generated_summary))) {
    warning("Reference/generated summaries do not share required comparison columns.")
    return(invisible(NULL))
  }
  gen <- generated_summary[, compare_cols, drop = FALSE]
  names(gen)[-(seq_along(key_cols))] <- paste0(names(gen)[-(seq_along(key_cols))], "_generated")
  ref_keep <- ref[, compare_cols, drop = FALSE]
  names(ref_keep)[-(seq_along(key_cols))] <- paste0(names(ref_keep)[-(seq_along(key_cols))], "_reference")
  cmp <- merge(ref_keep, gen, by = key_cols, all = TRUE, sort = FALSE)
  cmp$objective_min_seed_match <- cmp$objective_min_seed_reference == cmp$objective_min_seed_generated
  cmp$n_seeds_match <- cmp$n_seeds_reference == cmp$n_seeds_generated
  cmp$subcluster_k_match <- cmp$subcluster_k_reference == cmp$subcluster_k_generated
  write_csv(cmp, file.path(out_dir, "landscape_subcluster_reference_comparison.csv"))
  invisible(cmp)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  project_root <- normalizePath(path.expand(as_chr(argv$project_root, default_project_root())), mustWork = FALSE)
  out_dir_arg <- as_chr(argv$out_dir, "")
  if (!nzchar(out_dir_arg)) stop("--out_dir is required.", call. = FALSE)
  out_dir <- normalizePath(path.expand(out_dir_arg), mustWork = FALSE)
  invivo_run_dir <- normalizePath(path.expand(as_chr(argv$invivo_run_dir)), mustWork = TRUE)
  invitro_run_dir <- normalizePath(path.expand(as_chr(argv$invitro_run_dir)), mustWork = TRUE)
  invivo_parameter_table <- normalizePath(
    path.expand(as_chr(argv$invivo_parameter_table, default_invivo_parameter_table(project_root))),
    mustWork = FALSE
  )
  invitro_parameter_table <- normalizePath(
    path.expand(as_chr(argv$invitro_parameter_table %||% argv$parameter_table, default_invitro_parameter_table(project_root))),
    mustWork = FALSE
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  result_root <- normalizePath(path.expand(as_chr(argv$result_root, default_result_root(project_root))), mustWork = FALSE)
  analysis_root <- normalizePath(path.expand(as_chr(argv$analysis_root, landscape_analysis_root(out_dir))), mustWork = FALSE)
  seed_table_dir <- file.path(analysis_root, "SeedParameterTables")
  cluster_output_dir <- landscape_cluster_output_dir(analysis_root)
  reductions <- unique(vapply(as_char_vec(argv$reductions %||% argv$reduction, c("tsne", "umap")), normalize_reduction, character(1L)))
  max_seeds <- as_int(argv$max_seeds, NA_integer_)
  n_threads <- as_int(argv$n_threads, max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L)))
  umap_seed <- as_int(argv$umap_seed, 123L)
  tsne_seed <- as_int(argv$tsne_seed, 123L)
  n_neighbors <- as_int(argv$n_neighbors, 80L)
  min_dist <- as_num(argv$min_dist, 0.1)
  tsne_perplexity <- as_num(argv$tsne_perplexity, 30)
  tsne_theta <- as_num(argv$tsne_theta, 0.5)
  tsne_max_iter <- as_int(argv$tsne_max_iter, 1000L)
  cluster_seed <- as_int(argv$cluster_seed, 123L)
  cluster_k_min <- as_int(argv$cluster_k_min, 2L)
  cluster_k_max <- as_int(argv$cluster_k_max, 8L)
  silhouette_sample_n <- as_int(argv$cluster_silhouette_sample_n, 5000L)
  subcluster_seed <- as_int(argv$subcluster_seed, cluster_seed + 1000L)
  subcluster_k_min <- as_int(argv$subcluster_k_min, 2L)
  subcluster_k_max <- as_int(argv$subcluster_k_max, 6L)
  subcluster_min_n <- as_int(argv$subcluster_min_n, 6L)
  drop_invivo_initial <- as_bool(argv$drop_invivo_parameter_table_initial %||% argv$drop_parameter_table_initial, TRUE)
  drop_invitro_initial <- as_bool(argv$drop_invitro_parameter_table_initial %||% argv$drop_parameter_table_initial, TRUE)
  pairing_policy <- normalize_pairing_policy(as_chr(argv$pairing_policy, "cartesian_by_method"))
  deduplicate_pairs <- as_bool(argv$deduplicate_pairs, FALSE)
  reference_arg <- as_chr(argv$reference_subcluster_dir, "")
  reference_dir <- if (nzchar(reference_arg)) normalizePath(path.expand(reference_arg), mustWork = FALSE) else ""
  invivo_best_csv_arg <- as_chr(argv$invivo_best_csv, "")
  invivo_initial_csv_arg <- as_chr(argv$invivo_initial_csv, "")
  invitro_best_csv_arg <- as_chr(argv$invitro_best_csv, "")
  invitro_initial_csv_arg <- as_chr(argv$invitro_initial_csv, "")
  invivo_seed_space_tasks <- as_chr(argv$invivo_seed_space_tasks %||% argv$invivo_tasks_tsv, "")
  invitro_seed_space_tasks <- as_chr(argv$invitro_seed_space_tasks %||% argv$invitro_tasks_tsv, "")
  prepare_only <- as_bool(argv$prepare_only, FALSE)
  invivo_curve_filter <- as_bool(argv$invivo_curve_filter, FALSE)
  invivo_curve_class <- as_chr(argv$invivo_curve_class, "monotone_increasing")
  invivo_curve_class_table <- as_chr(argv$invivo_curve_class_table, "")

  message("Project root: ", project_root)
  message("Result root: ", result_root)
  message("Analysis root: ", analysis_root)
  message("Output root: ", out_dir)
  message("Reductions: ", paste(reductions, collapse = ", "))
  message("In vivo parameter table fallback: ", invivo_parameter_table)
  message("In vitro parameter table fallback: ", invitro_parameter_table)

  use_prebuilt_tables <- all(nzchar(c(invivo_best_csv_arg, invivo_initial_csv_arg, invitro_best_csv_arg, invitro_initial_csv_arg)))
  if (isTRUE(use_prebuilt_tables)) {
    invivo_tables <- list(
      best_csv = normalizePath(path.expand(invivo_best_csv_arg), mustWork = TRUE),
      initial_csv = normalizePath(path.expand(invivo_initial_csv_arg), mustWork = TRUE),
      seed_dirs = if (nzchar(invivo_seed_space_tasks)) {
        read_tsv(normalizePath(path.expand(invivo_seed_space_tasks), mustWork = TRUE))$seed_dir
      } else {
        list_seed_dirs(invivo_run_dir)
      }
    )
    invitro_tables <- list(
      best_csv = normalizePath(path.expand(invitro_best_csv_arg), mustWork = TRUE),
      initial_csv = normalizePath(path.expand(invitro_initial_csv_arg), mustWork = TRUE),
      seed_dirs = if (nzchar(invitro_seed_space_tasks)) {
        read_tsv(normalizePath(path.expand(invitro_seed_space_tasks), mustWork = TRUE))$seed_dir
      } else {
        list_seed_dirs(invitro_run_dir)
      }
    )
    message("Using prebuilt seed-space tables.")
    message("  invivo_best_csv: ", invivo_tables$best_csv)
    message("  invivo_initial_csv: ", invivo_tables$initial_csv)
    message("  invitro_best_csv: ", invitro_tables$best_csv)
    message("  invitro_initial_csv: ", invitro_tables$initial_csv)
  } else {
    invivo_tables <- build_seed_parameter_tables("invivo", invivo_run_dir, seed_table_dir, max_seeds = max_seeds, parameter_table_fallback = invivo_parameter_table)
    invitro_tables <- build_seed_parameter_tables("invitro", invitro_run_dir, seed_table_dir, max_seeds = max_seeds, parameter_table_fallback = invitro_parameter_table)
  }
  write_landscape_input_tables(
    out_dir = out_dir,
    invivo_tables = invivo_tables,
    invitro_tables = invitro_tables,
    invivo_run_dir = invivo_run_dir,
    invitro_run_dir = invitro_run_dir,
    analysis_root = analysis_root
  )
  pooled <- prepare_pooled_tables(
    invivo_best_csv = invivo_tables$best_csv,
    invivo_initial_csv = invivo_tables$initial_csv,
    invitro_best_csv = invitro_tables$best_csv,
    invitro_initial_csv = invitro_tables$initial_csv,
    drop_invivo_parameter_table_initial = drop_invivo_initial,
    drop_invitro_parameter_table_initial = drop_invitro_initial
  )
  feature_mat <- standardize_features(rbind(pooled$initial_features, pooled$best_features))
  write_csv(
    data.frame(parameter = colnames(feature_mat), preprocessing = "zscore", stringsAsFactors = FALSE),
    file.path(analysis_root, "pooled_feature_metadata.csv")
  )

  coordinate_paths <- list()
  for (reduction in reductions) {
    emb <- run_reduction_embedding(
      feature_mat,
      reduction = reduction,
      label = paste("pooled full", reduction),
      umap_seed = umap_seed,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      n_threads = n_threads,
      tsne_seed = tsne_seed,
      tsne_perplexity = tsne_perplexity,
      tsne_theta = tsne_theta,
      tsne_max_iter = tsne_max_iter
    )
    coordinate_paths[[reduction]] <- write_reduction_coordinate_table(reduction, emb, pooled, analysis_root)
  }

  feature_data <- load_pooled_best_feature_data(invivo_tables$best_csv, invitro_tables$best_csv)
  cluster_results <- lapply(reductions, function(reduction) {
    analyze_embedding(
      reduction = reduction,
      coordinate_csv = coordinate_paths[[reduction]],
      feature_data = feature_data,
      output_dir = cluster_output_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      silhouette_sample_n = silhouette_sample_n,
      subcluster_seed = subcluster_seed,
      subcluster_k_min = subcluster_k_min,
      subcluster_k_max = subcluster_k_max,
      subcluster_min_n = subcluster_min_n
    )
  })
  names(cluster_results) <- reductions
  all_seed_groups <- rbind_fill_plain(lapply(cluster_results, `[[`, "seed_groups"))
  all_subcluster_seed_groups <- rbind_fill_plain(lapply(cluster_results, `[[`, "subcluster_seed_groups"))
  all_subcluster_summary <- rbind_fill_plain(lapply(cluster_results, `[[`, "subcluster_summary"))
  combined_table_dir <- landscape_combined_table_dir(analysis_root)
  write_csv(all_seed_groups, file.path(combined_table_dir, "pooled_invivo_invitro_best_seed_groups_by_method.csv"))
  write_csv(all_subcluster_seed_groups, file.path(combined_table_dir, "pooled_invivo_invitro_best_subclusters_by_method.csv"))
  write_csv(all_subcluster_summary, file.path(combined_table_dir, "pooled_invivo_invitro_best_subcluster_summary_by_method.csv"))

  if (isTRUE(prepare_only)) {
    write_seed_plan_mode(
      out_dir,
      mode = "landscape_subcluster_prepared",
      warmup_pairs = 0L,
      pairing_policy = pairing_policy,
      deduplicate_pairs = deduplicate_pairs
    )
    write_validation_summary(reference_dir, all_subcluster_summary, out_dir)
    message("Prepared landscape subcluster outputs without finalizing warm-up pairs: ", combined_table_dir)
    return(invisible(list(seed_groups = all_seed_groups, subcluster_seed_groups = all_subcluster_seed_groups, subcluster_summary = all_subcluster_summary)))
  }

  finalize_landscape_pairs_core(
    out_dir = out_dir,
    project_root = project_root,
    invivo_run_dir = invivo_run_dir,
    invitro_run_dir = invitro_run_dir,
    analysis_root = analysis_root,
    result_root = result_root,
    reductions = reductions,
    umap_seed = umap_seed,
    tsne_seed = tsne_seed,
    cluster_seed = cluster_seed,
    subcluster_seed = subcluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    subcluster_k_min = subcluster_k_min,
    subcluster_k_max = subcluster_k_max,
    subcluster_min_n = subcluster_min_n,
    pairing_policy = pairing_policy,
    deduplicate_pairs = deduplicate_pairs,
    reference_dir = reference_dir,
    invivo_curve_filter = invivo_curve_filter,
    invivo_curve_class = invivo_curve_class,
    invivo_curve_class_table = invivo_curve_class_table,
    invitro_best_csv = invitro_tables$best_csv
  )
}

main_finalize_pairs <- function(argv) {
  project_root <- normalizePath(path.expand(as_chr(argv$project_root, default_project_root())), mustWork = FALSE)
  out_dir_arg <- as_chr(argv$out_dir, "")
  if (!nzchar(out_dir_arg)) stop("--out_dir is required.", call. = FALSE)
  out_dir <- normalizePath(path.expand(out_dir_arg), mustWork = FALSE)
  input_tables <- read_key_value_table(landscape_input_tables_path(out_dir))
  invivo_run_dir <- normalizePath(
    path.expand(as_chr(argv$invivo_run_dir, table_value(input_tables, "invivo_run_dir", ""))),
    mustWork = TRUE
  )
  invitro_run_dir <- normalizePath(
    path.expand(as_chr(argv$invitro_run_dir, table_value(input_tables, "invitro_run_dir", ""))),
    mustWork = TRUE
  )
  result_root <- normalizePath(path.expand(as_chr(argv$result_root, default_result_root(project_root))), mustWork = FALSE)
  analysis_root <- normalizePath(
    path.expand(as_chr(argv$analysis_root, table_value(input_tables, "analysis_root", landscape_analysis_root(out_dir)))),
    mustWork = FALSE
  )
  reductions <- unique(vapply(as_char_vec(argv$reductions %||% argv$reduction, c("tsne", "umap")), normalize_reduction, character(1L)))
  cluster_seed <- as_int(argv$cluster_seed, 123L)
  subcluster_seed <- as_int(argv$subcluster_seed, cluster_seed + 1000L)
  reference_arg <- as_chr(argv$reference_subcluster_dir, "")
  reference_dir <- if (nzchar(reference_arg)) normalizePath(path.expand(reference_arg), mustWork = FALSE) else ""
  finalize_landscape_pairs_core(
    out_dir = out_dir,
    project_root = project_root,
    invivo_run_dir = invivo_run_dir,
    invitro_run_dir = invitro_run_dir,
    analysis_root = analysis_root,
    result_root = result_root,
    reductions = reductions,
    umap_seed = as_int(argv$umap_seed, 123L),
    tsne_seed = as_int(argv$tsne_seed, 123L),
    cluster_seed = cluster_seed,
    subcluster_seed = subcluster_seed,
    cluster_k_min = as_int(argv$cluster_k_min, 2L),
    cluster_k_max = as_int(argv$cluster_k_max, 8L),
    subcluster_k_min = as_int(argv$subcluster_k_min, 2L),
    subcluster_k_max = as_int(argv$subcluster_k_max, 6L),
    subcluster_min_n = as_int(argv$subcluster_min_n, 6L),
    pairing_policy = normalize_pairing_policy(as_chr(argv$pairing_policy, "cartesian_by_method")),
    deduplicate_pairs = as_bool(argv$deduplicate_pairs, FALSE),
    reference_dir = reference_dir,
    invivo_curve_filter = as_bool(argv$invivo_curve_filter, TRUE),
    invivo_curve_class = as_chr(argv$invivo_curve_class, "monotone_increasing"),
    invivo_curve_class_table = as_chr(argv$invivo_curve_class_table, ""),
    invitro_best_csv = as_chr(argv$invitro_best_csv, table_value(input_tables, "invitro_best_csv", ""))
  )
}

main_build_seed_space_tasks <- function(argv) {
  project_root <- normalizePath(path.expand(as_chr(argv$project_root, default_project_root())), mustWork = FALSE)
  dataset <- match.arg(as_chr(argv$dataset), c("invivo", "invitro"))
  input_dir <- normalizePath(path.expand(as_chr(argv$input_dir %||% argv$run_dir)), mustWork = TRUE)
  out_dir <- normalizePath(path.expand(as_chr(argv$out_dir)), mustWork = FALSE)
  if (!nzchar(out_dir)) stop("--out_dir is required.", call. = FALSE)
  tables_dir <- normalizePath(
    path.expand(as_chr(argv$tables_dir, file.path(out_dir, "..", "landscape_subcluster", "SeedParameterTables"))),
    mustWork = FALSE
  )
  parameter_table_fallback <- as_chr(
    argv$parameter_table_fallback %||% argv$parameter_table,
    if (identical(dataset, "invivo")) default_invivo_parameter_table(project_root) else default_invitro_parameter_table(project_root)
  )
  build_seed_space_tasks(
    dataset = dataset,
    input_dir = input_dir,
    out_dir = out_dir,
    tables_dir = tables_dir,
    max_seeds = as_int(argv$max_seeds, NA_integer_),
    parameter_table_fallback = parameter_table_fallback
  )
}

main_run_seed_space_task <- function(argv) {
  tasks_tsv <- normalizePath(path.expand(as_chr(argv$tasks_tsv %||% argv$task_table)), mustWork = TRUE)
  task_id <- as_int(argv$task_id, as_int(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"), 1L))
  task_lookup_column <- as_chr(argv$task_lookup_column, "task_id")
  run_seed_space_task(
    tasks_tsv = tasks_tsv,
    task_id = task_id,
    task_lookup_column = task_lookup_column,
    skip_existing = as_bool(argv$skip_existing, TRUE)
  )
}

main_collect_seed_space_tables <- function(argv) {
  dataset <- match.arg(as_chr(argv$dataset), c("invivo", "invitro"))
  tasks_tsv <- normalizePath(path.expand(as_chr(argv$tasks_tsv %||% argv$task_table)), mustWork = TRUE)
  tables_dir <- normalizePath(path.expand(as_chr(argv$tables_dir)), mustWork = FALSE)
  if (!nzchar(tables_dir)) stop("--tables_dir is required.", call. = FALSE)
  collect_seed_space_tables(
    dataset = dataset,
    tasks_tsv = tasks_tsv,
    tables_dir = tables_dir
  )
}

dispatch_main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  stage <- tolower(gsub("-", "_", as_chr(argv$stage %||% argv$mode, "build_pairs"), fixed = TRUE))
  switch(
    stage,
    build_pairs = main(argv),
    pair = main(argv),
    prepare_landscape = {
      argv$prepare_only <- "TRUE"
      main(argv)
    },
    prepare = {
      argv$prepare_only <- "TRUE"
      main(argv)
    },
    finalize_pairs = main_finalize_pairs(argv),
    finalize = main_finalize_pairs(argv),
    build_seed_space_tasks = main_build_seed_space_tasks(argv),
    seed_space_tasks = main_build_seed_space_tasks(argv),
    run_seed_space_task = main_run_seed_space_task(argv),
    seed_space_task = main_run_seed_space_task(argv),
    collect_seed_space_tables = main_collect_seed_space_tables(argv),
    collect_seed_space = main_collect_seed_space_tables(argv),
    stop("Unknown --stage: ", stage, call. = FALSE)
  )
}

if (identical(environment(), globalenv())) {
  dispatch_main()
}
