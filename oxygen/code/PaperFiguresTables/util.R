#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (length(x) > 1L) return(x)
  if (is.na(x) || !nzchar(trimws(as.character(x)))) y else x
}

script_dir <- function() {
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
  if (length(frame_files) > 0L) return(dirname(frame_files[[length(frame_files)]]))
  normalizePath(getwd(), mustWork = FALSE)
}

default_oxygen_dir <- function() {
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cwd_oxygen <- file.path(cwd, "oxygen")
  if (dir.exists(file.path(cwd_oxygen, "code", "PaperFiguresTables"))) {
    return(cwd_oxygen)
  }

  sd <- script_dir()
  if (identical(basename(sd), "PaperFiguresTables") && identical(basename(dirname(sd)), "code")) {
    return(dirname(dirname(sd)))
  }
  dirname(dirname(sd))
}

default_paperfigures_dir <- function() {
  file.path(default_oxygen_dir(), "results", "PaperFiguresTables")
}

DEFAULT_INVIVO_INPUT_DIR <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed"
DEFAULT_INVITRO_INPUT_DIR <- "/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed"
DEFAULT_OUTPUT_DIR <- default_paperfigures_dir()

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  y <- tolower(trimws(as.character(x[[1]])))
  if (!nzchar(y)) return(default)
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

as_int <- function(x, default = NA_integer_) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  y <- suppressWarnings(as.integer(x[[1]]))
  if (!is.finite(y) || is.na(y)) default else y
}

as_num <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  y <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(y) || is.na(y)) default else y
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
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

metric_map <- function(path) {
  tab <- read_tsv(path)
  if (!all(c("metric", "value") %in% names(tab))) {
    stop("Summary file is not a metric/value table: ", path)
  }
  stats::setNames(as.character(tab$value), as.character(tab$metric))
}

metric_value <- function(metrics, keys, default = NA_character_) {
  hit <- keys[keys %in% names(metrics)]
  if (!length(hit)) return(default)
  val <- metrics[[hit[[1]]]]
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
  dirs[order(seed_from_dir(dirs))]
}

read_best_params <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  if (!file.exists(path)) stop("Missing best_params.tsv: ", path)
  tab <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv must contain parameter/value columns: ", path)
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

transform_param_slot <- function(value, transform, param_symbol, slot_label) {
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || is.na(value)) stop("Non-finite ", slot_label, " for parameter ", param_symbol)
  if (identical(transform, "log10")) {
    if (value <= 0) stop(param_symbol, " ", slot_label, " must be > 0 for log10 transform.")
    return(log10(value))
  }
  value
}

convert_invitro_parameter_table <- function(tab, path) {
  req <- c("param_symbol", "init_value", "lower_bound", "upper_bound")
  missing <- setdiff(req, names(tab))
  if (length(missing)) stop("parameter_table_input.csv missing columns: ", paste(missing, collapse = ", "), " in ", path)
  estimate_col <- if ("use_invitro_fit" %in% names(tab)) "use_invitro_fit" else "estimate"
  if (!estimate_col %in% names(tab)) {
    stop("parameter_table_input.csv must contain use_invitro_fit or estimate in ", path)
  }

  tab$param_symbol <- trimws(as.character(tab$param_symbol))
  map <- invitro_parameter_transform_map()
  tab <- merge(map, tab, by = "param_symbol", all.x = FALSE, all.y = FALSE, sort = FALSE)
  tab <- tab[match(map$param_symbol, tab$param_symbol, nomatch = 0L), , drop = FALSE]
  if (!nrow(tab)) stop("No supported in vitro optimizer parameters found in ", path)

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

read_param_table <- function(seed_dir, dataset = "invivo") {
  candidates <- c(
    file.path(seed_dir, "parameter_table.csv"),
    file.path(dirname(seed_dir), "parameter_table.csv"),
    file.path(seed_dir, "parameter_table_input.csv"),
    file.path(dirname(seed_dir), "parameter_table_input.csv"),
    file.path(seed_dir, "parameter_table_input_invitro.csv"),
    file.path(dirname(seed_dir), "parameter_table_input_invitro.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) stop("Missing parameter_table.csv for ", seed_dir)
  path <- existing[[1L]]
  tab <- read_csv_plain(path)
  req <- c("param_name", "estimate", "init_value", "lower_bound", "upper_bound", "param_prototype")
  missing <- setdiff(req, names(tab))
  if (length(missing)) {
    if (identical(dataset, "invitro") && "param_symbol" %in% names(tab)) {
      return(convert_invitro_parameter_table(tab, path))
    }
    stop("parameter_table.csv missing columns: ", paste(missing, collapse = ", "), " in ", path)
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
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  u <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  span <- as.numeric(upper - lower)
  out <- sweep(u, 2L, span, `*`)
  out <- sweep(out, 2L, as.numeric(lower), `+`)
  colnames(out) <- names(lower)
  out
}

sample_truncnorm_box <- function(n, center, lower, upper, sigma_frac = 0.1) {
  n <- as.integer(n)
  d <- length(lower)
  if (is.na(n) || n <= 0L || d <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = d))
  }
  sigma_frac <- as.numeric(sigma_frac)
  if (!is.finite(sigma_frac) || sigma_frac <= 0) sigma_frac <- 0.1
  sd_vec <- pmax((as.numeric(upper) - as.numeric(lower)) * sigma_frac, 1e-12)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  out <- sweep(z, 2L, sd_vec, `*`)
  out <- sweep(out, 2L, as.numeric(center), `+`)
  out <- sweep(out, 2L, as.numeric(lower), pmax)
  out <- sweep(out, 2L, as.numeric(upper), pmin)
  colnames(out) <- names(lower)
  out
}

build_de_initialpop <- function(np,
                                lower,
                                upper,
                                init_use,
                                mode = "uniform",
                                uniform_frac = 0.3,
                                sigma_frac = 0.1) {
  np <- as.integer(np)
  if (is.na(np) || np < 1L) stop("NP must be >= 1")
  mode <- tolower(trimws(as.character(mode)))
  if (!mode %in% c("uniform", "hybrid")) stop("de_init_mode must be uniform or hybrid")

  pop <- sample_uniform_box(np, lower, upper)
  init_vec <- as.numeric(init_use[names(lower)])
  names(init_vec) <- names(lower)
  if (any(!is.finite(init_vec))) stop("Initial vector has missing/non-finite values")
  init_vec <- pmin(pmax(init_vec, lower), upper)
  pop[1L, ] <- init_vec

  rem <- as.integer(np - 1L)
  if (rem <= 0L) return(pop)

  if (identical(mode, "uniform")) {
    pop[2L:np, ] <- sample_uniform_box(rem, lower, upper)
    return(pop)
  }

  uniform_frac <- as.numeric(uniform_frac)
  if (!is.finite(uniform_frac) || uniform_frac < 0 || uniform_frac > 1) {
    stop("de_init_uniform_frac must be in [0,1]")
  }
  n_uniform <- as.integer(round(rem * uniform_frac))
  n_uniform <- max(0L, min(rem, n_uniform))
  n_local <- as.integer(rem - n_uniform)
  cursor <- 2L
  if (n_local > 0L) {
    idx <- cursor:(cursor + n_local - 1L)
    pop[idx, ] <- sample_truncnorm_box(n_local, init_vec, lower, upper, sigma_frac)
    cursor <- cursor + n_local
  }
  if (n_uniform > 0L) {
    idx <- cursor:(cursor + n_uniform - 1L)
    pop[idx, ] <- sample_uniform_box(n_uniform, lower, upper)
  }
  if (rem > 1L) {
    ord <- sample.int(rem, size = rem, replace = FALSE)
    pop[2L:np, ] <- pop[1L + ord, , drop = FALSE]
  }
  pop
}

inverse_transform_column <- function(param_name, x) {
  if (startsWith(param_name, "log10_")) {
    return(10^as.numeric(x))
  }
  as.numeric(x)
}

initial_population_natural <- function(seed_dir, seed, metrics, best_names, dataset = "invivo") {
  param_table <- read_param_table(seed_dir, dataset = dataset)
  opt_tab <- param_table[param_table$estimate, , drop = FALSE]
  if (!nrow(opt_tab)) stop("No estimated parameters in parameter table for ", seed_dir)

  param_names <- as.character(opt_tab$param_name)
  init <- stats::setNames(as.numeric(opt_tab$init_value), param_names)
  lower <- stats::setNames(as.numeric(opt_tab$lower_bound), param_names)
  upper <- stats::setNames(as.numeric(opt_tab$upper_bound), param_names)
  np <- read_np(metrics, param_table)

  set.seed(as.integer(seed))
  pop_t <- if (identical(dataset, "invitro")) {
    sample_uniform_box(np, lower, upper)
  } else {
    build_de_initialpop(
      np = np,
      lower = lower,
      upper = upper,
      init_use = init,
      mode = "uniform"
    )
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
    if (nm %in% names(best_vals)) {
      out[[nm]] <- as.numeric(best_vals[[nm]])
    } else {
      out[[nm]] <- NA_real_
    }
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
  for (nm in intersect(best_names, names(vals))) {
    out[[nm]] <- as.numeric(vals[[nm]])
  }
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
  } else if (identical(dataset, "invitro")) {
    out$objective_growth <- negative_metric_numeric(metrics, "growth_loglik", seed, "objective_growth")
    out$objective_ploidy <- negative_metric_numeric(metrics, "ploidy_loglik", seed, "objective_ploidy")
    out$objective_flow <- negative_metric_numeric(metrics, "flow_loglik", seed, "objective_flow")
  }
  iter_completed <- suppressWarnings(as.integer(metric_value(
    metrics,
    c("optimizer_iter_completed", "deoptim_iter_completed", "deoptim_stop_iteration", "iter_completed"),
    default = NA_character_
  )))
  if (!is.finite(iter_completed) || is.na(iter_completed)) {
    iter_completed <- suppressWarnings(as.integer(metric_value(metrics, "itermax", default = NA_character_)))
    warning("Using itermax fallback for DEoptim completed iteration in seed", seed)
  }
  out$deoptim_iter_completed <- iter_completed
  if (identical(dataset, "invivo")) {
    out$pred1000_both_gt44 <- read_pred1000_both_gt44(seed_dir, seed)
    out <- cbind(out, read_pred1000_ploidy_ratios(seed_dir, seed))
  }
  out$seed <- as.integer(seed)
  out
}

seed_prediction_path <- function(seed_dir, filename) {
  candidates <- c(
    file.path(seed_dir, "viz", filename),
    file.path(seed_dir, "viz", "invivo", filename)
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1L]] else candidates[[1L]]
}

prediction_value_column <- function(tab) {
  cols <- c("weighted_mean_N", "weighted_mean_endpoint", "weighted_mean_ploidy")
  hit <- cols[cols %in% names(tab)]
  if (length(hit)) hit[[1L]] else NA_character_
}

prediction_threshold <- function(tab, value_col) {
  if (identical(value_col, "weighted_mean_N")) return(44)
  vals <- suppressWarnings(as.numeric(tab[[value_col]]))
  if (any(is.finite(vals) & vals > 10, na.rm = TRUE)) return(44)
  2
}

prediction_initial_ploidy <- function(tab, value_col, cohort) {
  cohort <- as.character(cohort)
  values <- suppressWarnings(as.numeric(tab[[value_col]]))
  start_with <- if ("start_with" %in% names(tab)) tolower(as.character(tab$start_with)) else character(0)
  chromosome_scale <- identical(value_col, "weighted_mean_N") ||
    any(grepl("chr", start_with), na.rm = TRUE) ||
    any(is.finite(values) & values > 10, na.rm = TRUE)
  if (isTRUE(chromosome_scale)) {
    if (identical(cohort, "4N")) 88 else 44
  } else {
    if (identical(cohort, "4N")) 4 else 2
  }
}

read_pred1000_ploidy_ratios <- function(seed_dir, seed, target_day = 1000) {
  out <- data.frame(
    pred1000_ploidy_ratio_2N = NA_real_,
    pred1000_ploidy_ratio_4N = NA_real_,
    stringsAsFactors = FALSE
  )
  path <- seed_prediction_path(seed_dir, "predict_ploidy_weighted_mean_0_1000day.tsv")
  if (!file.exists(path)) {
    warning("Missing day1000 ploidy prediction file for seed", seed, ": ", path)
    return(out)
  }

  tab <- tryCatch(read_tsv(path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab)) {
    warning("Could not read day1000 ploidy prediction file for seed", seed)
    return(out)
  }

  required <- c("cohort", "day")
  if (!all(required %in% names(tab))) {
    warning("Day1000 ploidy prediction file missing cohort/day columns for seed", seed)
    return(out)
  }
  value_col <- prediction_value_column(tab)
  if (is.na(value_col)) {
    warning("Day1000 ploidy prediction file missing weighted mean columns for seed", seed)
    return(out)
  }

  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab[[value_col]] <- suppressWarnings(as.numeric(tab[[value_col]]))
  target_day <- suppressWarnings(as.numeric(target_day))
  rows <- tab[
    is.finite(tab$day) &
      abs(tab$day - target_day) <= 1e-8 &
      as.character(tab$cohort) %in% c("2N", "4N") &
      is.finite(tab[[value_col]]),
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) {
    warning("No valid day1000 2N/4N ploidy prediction rows for seed", seed)
    return(out)
  }

  for (cohort_id in c("2N", "4N")) {
    cohort_rows <- rows[as.character(rows$cohort) == cohort_id, , drop = FALSE]
    if (!nrow(cohort_rows)) next
    denom <- prediction_initial_ploidy(cohort_rows, value_col, cohort_id)
    if (!is.finite(denom) || denom <= 0) next
    col <- paste0("pred1000_ploidy_ratio_", cohort_id)
    out[[col]] <- mean(cohort_rows[[value_col]], na.rm = TRUE) / denom
  }
  out
}

read_pred1000_both_gt44 <- function(seed_dir, seed, target_day = 1000) {
  path <- seed_prediction_path(seed_dir, "predict_ploidy_weighted_mean_0_1000day.tsv")
  if (!file.exists(path)) {
    warning("Missing day1000 ploidy prediction file for seed", seed, ": ", path)
    return(FALSE)
  }

  tab <- tryCatch(read_tsv(path), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab)) {
    warning("Could not read day1000 ploidy prediction file for seed", seed)
    return(FALSE)
  }

  required <- c("cohort", "day")
  if (!all(required %in% names(tab))) {
    warning("Day1000 ploidy prediction file missing cohort/day columns for seed", seed)
    return(FALSE)
  }
  value_col <- prediction_value_column(tab)
  if (is.na(value_col)) {
    warning("Day1000 ploidy prediction file missing weighted mean columns for seed", seed)
    return(FALSE)
  }

  tab$day <- suppressWarnings(as.numeric(tab$day))
  tab[[value_col]] <- suppressWarnings(as.numeric(tab[[value_col]]))
  target_day <- suppressWarnings(as.numeric(target_day))
  rows <- tab[
    is.finite(tab$day) &
      abs(tab$day - target_day) <= 1e-8 &
      as.character(tab$cohort) %in% c("2N", "4N") &
      is.finite(tab[[value_col]]),
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) {
    warning("No valid day1000 2N/4N ploidy prediction rows for seed", seed)
    return(FALSE)
  }

  threshold <- prediction_threshold(rows, value_col)
  values <- split(rows[[value_col]], as.character(rows$cohort))
  isTRUE(
    all(c("2N", "4N") %in% names(values)) &&
      all(values[["2N"]] > threshold) &&
      all(values[["4N"]] > threshold)
  )
}

write_csv <- function(df, path) {
  utils::write.csv(df, file = path, quote = FALSE, row.names = FALSE)
  message("Wrote ", path, " [", nrow(df), " x ", ncol(df), "]")
}

normalize_dataset <- function(dataset) {
  dataset <- tolower(trimws(as.character(dataset %||% "invivo")))
  if (!dataset %in% c("invivo", "invitro")) stop("dataset must be invivo or invitro.")
  dataset
}

dataset_label <- function(dataset) {
  dataset <- normalize_dataset(dataset)
  if (identical(dataset, "invivo")) "in vivo" else "in vitro"
}

default_paper_figures_tables_dir <- function() {
  file.path(default_oxygen_dir(), "results", "PaperFiguresTables")
}

paper_dataset_dir <- function(dataset, root_dir = default_paper_figures_tables_dir()) {
  file.path(normalizePath(path.expand(root_dir), mustWork = FALSE), normalize_dataset(dataset))
}

paper_tables_dir <- function(dataset, root_dir = default_paper_figures_tables_dir()) {
  file.path(paper_dataset_dir(dataset, root_dir = root_dir), "Tables")
}

paper_figures_dir <- function(dataset, root_dir = default_paper_figures_tables_dir()) {
  file.path(paper_dataset_dir(dataset, root_dir = root_dir), "Figures")
}

paper_tables_wclusters_dir <- function(dataset, root_dir = default_paper_figures_tables_dir()) {
  file.path(paper_dataset_dir(dataset, root_dir = root_dir), "TablesWclusters")
}

paper_figures_wclusters_dir <- function(dataset, root_dir = default_paper_figures_tables_dir()) {
  file.path(paper_dataset_dir(dataset, root_dir = root_dir), "FiguresWclusters")
}

default_dataset_input_dir <- function(dataset) {
  dataset <- normalize_dataset(dataset)
  if (identical(dataset, "invitro")) DEFAULT_INVITRO_INPUT_DIR else DEFAULT_INVIVO_INPUT_DIR
}

paper_generate_umap_tables <- function(dataset = "invivo",
                                       input_dir = default_dataset_input_dir(dataset),
                                       tables_dir = paper_tables_dir(dataset),
                                       max_seeds = NA_integer_,
                                       write_best = TRUE,
                                       write_initial = TRUE) {
  dataset <- normalize_dataset(dataset)
  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  max_seeds <- as_int(max_seeds, NA_integer_)

  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  if (!write_best && !write_initial) stop("Nothing to write: set write_best and/or write_initial to TRUE.")

  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No seed directories found under: ", input_dir)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
  }

  best_vectors <- list()
  for (seed_dir in seed_dirs) {
    seed <- seed_from_dir(seed_dir)
    best_vectors[[as.character(seed)]] <- read_best_params(seed_dir)
  }
  best_names <- names(best_vectors[[1L]])
  for (vals in best_vectors[-1L]) {
    best_names <- c(best_names, setdiff(names(vals), best_names))
  }

  message("Found ", length(seed_dirs), " ", dataset_label(dataset), " seed directories.")
  message("Using ", length(best_names), " natural parameter columns.")

  best_rows <- if (write_best) vector("list", length(seed_dirs)) else NULL
  init_rows <- if (write_initial) vector("list", length(seed_dirs)) else NULL
  for (i in seq_along(seed_dirs)) {
    seed_dir <- seed_dirs[[i]]
    seed <- seed_from_dir(seed_dir)
    summary_path <- file.path(seed_dir, "fit_summary.tsv")
    if (!file.exists(summary_path)) stop("Missing fit_summary.tsv: ", summary_path)
    metrics <- metric_map(summary_path)

    if (write_best) {
      best_rows[[i]] <- best_params_row(seed_dir, seed, best_names, dataset = dataset)
    }
    if (write_initial) {
      init_rows[[i]] <- initial_population_natural(seed_dir, seed, metrics, best_names, dataset = dataset)
    }

    if (i %% 25L == 0L || i == length(seed_dirs)) {
      message("Processed ", i, "/", length(seed_dirs), " seeds.")
    }
  }

  out <- list()
  if (write_initial) {
    init_df <- do.call(rbind, init_rows)
    init_df <- init_df[, c(best_names, "seed"), drop = FALSE]
    out$initial_population <- file.path(tables_dir, paste0(dataset, "_deoptim_initial_population.csv"))
    write_csv(init_df, out$initial_population)
  }
  if (write_best) {
    best_df <- do.call(rbind, best_rows)
    extra_cols <- if (identical(dataset, "invivo")) {
      c(
        "objective", "objective_ploidy", "objective_burden", "deoptim_iter_completed",
        "pred1000_both_gt44", "pred1000_ploidy_ratio_2N", "pred1000_ploidy_ratio_4N",
        "seed"
      )
    } else {
      c("objective", "objective_growth", "objective_ploidy", "objective_flow", "deoptim_iter_completed", "seed")
    }
    best_df <- best_df[, c(best_names, extra_cols), drop = FALSE]
    out$best_params <- file.path(tables_dir, paste0(dataset, "_best_params_by_seed.csv"))
    write_csv(best_df, out$best_params)
  }

  invisible(out)
}

read_objective_from_summary <- function(seed_dir) {
  path <- file.path(seed_dir, "fit_summary.tsv")
  if (!file.exists(path)) return(NA_real_)
  metrics <- metric_map(path)
  metric_numeric(
    metrics,
    c("objective", "objective_total", "optimizer_local_objective", "optimizer_deoptim_objective"),
    seed = seed_from_dir(seed_dir),
    label = "objective"
  )
}

attach_objective <- function(best_df, objective_seed_dir) {
  if ("objective" %in% names(best_df)) {
    best_df$objective <- suppressWarnings(as.numeric(best_df$objective))
    return(best_df)
  }
  if (is.null(objective_seed_dir) || !dir.exists(objective_seed_dir)) {
    stop(
      "Best-parameter CSV has no objective column, and objective_seed_dir does not exist: ",
      objective_seed_dir,
      "\nProvide --objective_seed_dir=/path/to/fit_*_O2_buffering_500seed or add objective to the best-parameter CSV."
    )
  }
  if (!"seed" %in% names(best_df)) stop("Best-parameter CSV must contain a seed column.")
  best_df$objective <- vapply(
    best_df$seed,
    function(seed) read_objective_from_summary(file.path(objective_seed_dir, paste0("seed", as.integer(seed)))),
    numeric(1)
  )
  if (any(!is.finite(best_df$objective))) {
    bad <- best_df$seed[!is.finite(best_df$objective)]
    stop("Missing/non-finite objective values for seeds: ", paste(head(bad, 20), collapse = ", "))
  }
  best_df
}

coerce_logical_column <- function(x, label) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (any(is.na(out))) {
    stop("Column ", label, " contains values that cannot be converted to TRUE/FALSE.")
  }
  out
}

umap_parameter_set <- function(dataset) {
  dataset <- normalize_dataset(dataset)
  if (identical(dataset, "invitro")) {
    return(c(
      "O2_crit",
      "mu_hp",
      "p_misseg",
      "k_o_mis",
      "buffer_smax",
      "buffer_beta",
      "buffer_n_exp",
      "n_O",
      "alpha_o2",
      "gamma_growth",
      "lam_max",
      "p_mis_base",
      "p_wgd",
      "gamma_mu"
    ))
  }
  c(
    "lam_max",
    "p_mis_base",
    "p_misseg",
    "k_o_mis",
    "buffer_smax",
    "buffer_beta",
    "buffer_n_exp",
    "p_wgd",
    "o2_S0",
    "kappa_O",
    "eta_o2",
    "alpha_o2",
    "gamma_growth",
    "mu_hp",
    "gamma_mu",
    "O2_crit",
    "n_O",
    "k_clear"
  )
}

umap_log10_parameter_set <- function(dataset) {
  intersect(
    c(
      "lam_max",
      "p_mis_base",
      "p_misseg",
      "k_o_mis",
      "buffer_beta",
      "buffer_n_exp",
      "p_wgd",
      "o2_S0",
      "kappa_O",
      "eta_o2",
      "alpha_o2",
      "mu_hp",
      "O2_crit",
      "k_clear"
    ),
    umap_parameter_set(dataset)
  )
}

transform_umap_features <- function(df, params, log10_params) {
  missing <- setdiff(params, names(df))
  if (length(missing)) stop("Input table is missing UMAP parameter columns: ", paste(missing, collapse = ", "))

  out <- as.data.frame(lapply(params, function(param) {
    vals <- suppressWarnings(as.numeric(df[[param]]))
    if (param %in% log10_params) {
      if (any(!is.finite(vals) | vals <= 0)) {
        stop("Cannot log10-transform non-positive/non-finite values for parameter: ", param)
      }
      vals <- log10(vals)
    }
    vals
  }), check.names = FALSE)
  names(out) <- params

  if (any(!is.finite(as.matrix(out)))) stop("UMAP feature matrix contains non-finite values.")
  out
}

standardize_features <- function(x) {
  scaled <- scale(as.matrix(x))
  zero_sd <- !is.finite(attr(scaled, "scaled:scale")) | attr(scaled, "scaled:scale") == 0
  if (any(zero_sd)) {
    stop("UMAP feature columns have zero/non-finite SD after transformation: ", paste(names(x)[zero_sd], collapse = ", "))
  }
  scaled
}

drop_parameter_table_initial_rows <- function(initial_df) {
  row_in_seed <- ave(seq_len(nrow(initial_df)), initial_df$seed, FUN = seq_along)
  keep <- row_in_seed != 1L
  dropped <- sum(!keep)
  message("Dropped ", dropped, " parameter-table initial-value points from initial population.")
  initial_df[keep, , drop = FALSE]
}

default_pca_output_prefix <- function(output_prefix) {
  if (grepl("_umap$", output_prefix)) {
    return(sub("_umap$", "_pca_umap", output_prefix))
  }
  paste0(output_prefix, "_pca_umap")
}

default_sampled_output_prefix <- function(output_prefix, n_sample) {
  suffix <- paste0("_sampled", as.integer(n_sample), "_umap")
  if (grepl("_umap$", output_prefix)) {
    return(sub("_umap$", suffix, output_prefix))
  }
  paste0(output_prefix, suffix)
}

run_umap_embedding <- function(feature_mat,
                               label,
                               umap_seed,
                               n_neighbors,
                               min_dist,
                               n_threads) {
  message("Running ", label, " UMAP with ", nrow(feature_mat), " points and ", ncol(feature_mat), " input dimensions.")
  set.seed(umap_seed)
  emb <- uwot::umap(
    feature_mat,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "euclidean",
    n_threads = n_threads,
    n_sgd_threads = 1L,
    ret_model = FALSE,
    verbose = TRUE
  )
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb
}

build_plot_data <- function(emb, initial_df, best_df, reduction_label, shape_by_pred = TRUE) {
  n_initial <- nrow(initial_df)
  n_best <- nrow(best_df)
  plot_data <- data.frame(
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    point_type = c(rep("initial", n_initial), rep("best", n_best)),
    seed = c(as.integer(initial_df$seed), as.integer(best_df$seed)),
    objective = c(rep(NA_real_, n_initial), best_df$objective),
    reduction = reduction_label
  )
  if (isTRUE(shape_by_pred)) {
    if (!"pred1000_both_gt44" %in% names(best_df)) {
      stop("Best-parameter CSV must contain pred1000_both_gt44 when shape_by_pred=TRUE.")
    }
    plot_data$pred1000_both_gt44 <- c(rep(NA, n_initial), best_df$pred1000_both_gt44)
    plot_data$pred1000_both_gt44_label <- factor(
      ifelse(is.na(plot_data$pred1000_both_gt44), NA_character_, as.character(plot_data$pred1000_both_gt44)),
      levels = c("FALSE", "TRUE")
    )
  } else {
    plot_data$pred1000_both_gt44 <- NA
    plot_data$pred1000_both_gt44_label <- factor(NA_character_, levels = c("FALSE", "TRUE"))
  }
  plot_data
}

square_umap_limits <- function(plot_data, pad_frac = 0.05) {
  x <- suppressWarnings(as.numeric(plot_data$UMAP1))
  y <- suppressWarnings(as.numeric(plot_data$UMAP2))
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) stop("UMAP plot data must contain finite UMAP1/UMAP2 coordinates.")

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

build_umap_plot <- function(plot_data, initial_size = 0.22, best_size = 1.2, shape_by_pred = TRUE) {
  lims <- square_umap_limits(plot_data)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_data[plot_data$point_type == "initial", , drop = FALSE],
      ggplot2::aes(x = UMAP1, y = UMAP2),
      color = "grey48",
      alpha = 0.34,
      size = initial_size,
      stroke = 0
    )

  if (isTRUE(shape_by_pred)) {
    p <- p +
      ggplot2::geom_point(
        data = plot_data[plot_data$point_type == "best", , drop = FALSE],
        ggplot2::aes(x = UMAP1, y = UMAP2, color = objective, shape = pred1000_both_gt44_label),
        alpha = 0.95,
        size = best_size,
        stroke = 0
      ) +
      ggplot2::scale_shape_manual(
        name = "Pred1000 > 44",
        values = c("FALSE" = 16, "TRUE" = 18),
        labels = c("FALSE", "TRUE")
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        data = plot_data[plot_data$point_type == "best", , drop = FALSE],
        ggplot2::aes(x = UMAP1, y = UMAP2, color = objective),
        shape = 16,
        alpha = 0.95,
        size = best_size,
        stroke = 0
      )
  }

  p +
    ggplot2::scale_color_gradient(
      name = "Objective",
      low = "#2C7BB6",
      high = "#FDE725",
      guide = ggplot2::guide_colorbar(barheight = ggplot2::unit(35, "mm"))
    ) +
    ggplot2::coord_equal(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

cluster_output_subdirs <- function() {
  c(
    umap_coordinates = "ByUMAPCoordinates",
    input_features = "ByInputFeatures"
  )
}

relabel_clusters_by_umap <- function(cluster_num, plot_data) {
  centers <- stats::aggregate(
    cbind(UMAP1, UMAP2) ~ cluster,
    data = data.frame(
      cluster = as.integer(cluster_num),
      UMAP1 = suppressWarnings(as.numeric(plot_data$UMAP1)),
      UMAP2 = suppressWarnings(as.numeric(plot_data$UMAP2))
    ),
    FUN = stats::median
  )
  centers <- centers[order(centers$UMAP1, centers$UMAP2), , drop = FALSE]
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
                                   seed = 123L,
                                   k_min = 2L,
                                   k_max = 8L,
                                   silhouette_sample_n = 5000L,
                                   nstart = 10L,
                                   iter.max = 100L) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Required R package is not installed: cluster")
  }
  mat <- as.matrix(basis_mat)
  storage.mode(mat) <- "double"
  if (nrow(mat) != nrow(plot_data)) {
    stop("Cluster feature matrix row count does not match UMAP plot data for ", cluster_source, ".")
  }
  if (!nrow(mat) || !ncol(mat)) stop("Cluster feature matrix is empty for ", cluster_source, ".")
  if (any(!is.finite(mat))) stop("Cluster feature matrix contains non-finite values for ", cluster_source, ".")

  n <- nrow(mat)
  silhouette_sample_n <- as.integer(silhouette_sample_n)
  if (!is.finite(silhouette_sample_n) || is.na(silhouette_sample_n) || silhouette_sample_n < 2L) {
    silhouette_sample_n <- min(5000L, n)
  }
  sample_n <- min(n, silhouette_sample_n)
  if (n < 3L || sample_n < 3L) {
    return(single_cluster_assignment(n, cluster_source, sample_n))
  }

  k_min <- max(2L, as.integer(k_min))
  k_max <- min(as.integer(k_max), sample_n - 1L)
  set.seed(as.integer(seed))
  sample_idx <- if (n > sample_n) sort(sample.int(n, sample_n)) else seq_len(n)
  sample_mat <- mat[sample_idx, , drop = FALSE]
  unique_sample_n <- nrow(unique(sample_mat))
  k_max <- min(k_max, unique_sample_n - 1L)
  if (!is.finite(k_max) || is.na(k_max) || k_max < k_min) {
    return(single_cluster_assignment(n, cluster_source, sample_n))
  }

  k_values <- seq.int(k_min, k_max)
  sample_dist <- stats::dist(sample_mat)
  scores <- rep(NA_real_, length(k_values))
  centers_by_k <- vector("list", length(k_values))

  for (i in seq_along(k_values)) {
    k <- k_values[[i]]
    set.seed(as.integer(seed) + k)
    km <- try(
      stats::kmeans(
        sample_mat,
        centers = k,
        nstart = as.integer(nstart),
        iter.max = as.integer(iter.max)
      ),
      silent = TRUE
    )
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
  if (!any(valid)) {
    return(single_cluster_assignment(n, cluster_source, sample_n))
  }

  selected_i <- which(valid)[which.max(summary$average_silhouette[valid])]
  selected_k <- summary$k[[selected_i]]
  summary$selected[[selected_i]] <- TRUE

  final_km <- try(
    stats::kmeans(
      mat,
      centers = centers_by_k[[selected_i]],
      iter.max = as.integer(iter.max),
      algorithm = "Lloyd"
    ),
    silent = TRUE
  )
  if (inherits(final_km, "try-error")) {
    set.seed(as.integer(seed) + selected_k + 1000L)
    final_km <- stats::kmeans(
      mat,
      centers = selected_k,
      nstart = as.integer(nstart),
      iter.max = as.integer(iter.max)
    )
  }

  cluster_num <- relabel_clusters_by_umap(final_km$cluster, plot_data)
  list(
    cluster_num = cluster_num,
    cluster_id = sprintf("C%02d", cluster_num),
    summary = summary
  )
}

cluster_palette <- function(cluster_ids) {
  cluster_ids <- sort(unique(as.character(cluster_ids)))
  stats::setNames(grDevices::hcl.colors(length(cluster_ids), palette = "Dark 3"), cluster_ids)
}

cluster_hull_data <- function(plot_data, cluster_col = "cluster_id", expand = 1.035) {
  if (!cluster_col %in% names(plot_data)) return(data.frame())
  cluster_ids <- sort(unique(as.character(plot_data[[cluster_col]])))
  cluster_ids <- cluster_ids[nzchar(cluster_ids) & !is.na(cluster_ids)]
  if (!length(cluster_ids)) return(data.frame())

  all_span <- max(diff(range(plot_data$UMAP1, finite = TRUE)), diff(range(plot_data$UMAP2, finite = TRUE)))
  point_radius <- if (is.finite(all_span) && all_span > 0) all_span * 0.012 else 0.1

  hulls <- lapply(cluster_ids, function(cluster_id) {
    d <- plot_data[as.character(plot_data[[cluster_col]]) == cluster_id, c("UMAP1", "UMAP2"), drop = FALSE]
    d <- d[is.finite(d$UMAP1) & is.finite(d$UMAP2), , drop = FALSE]
    d <- unique(d)
    if (nrow(d) >= 3L) {
      idx <- grDevices::chull(d$UMAP1, d$UMAP2)
      h <- d[c(idx, idx[[1]]), , drop = FALSE]
      center <- colMeans(h[, c("UMAP1", "UMAP2"), drop = FALSE])
      h$UMAP1 <- center[[1]] + (h$UMAP1 - center[[1]]) * expand
      h$UMAP2 <- center[[2]] + (h$UMAP2 - center[[2]]) * expand
    } else if (nrow(d) >= 1L) {
      theta <- seq(0, 2 * pi, length.out = 33L)
      h <- data.frame(
        UMAP1 = d$UMAP1[[1]] + point_radius * cos(theta),
        UMAP2 = d$UMAP2[[1]] + point_radius * sin(theta)
      )
    } else {
      return(NULL)
    }
    h$cluster_id <- cluster_id
    h
  })
  out <- do.call(rbind, Filter(Negate(is.null), hulls))
  row.names(out) <- NULL
  out
}

cluster_label_data <- function(plot_data, cluster_col = "cluster_id") {
  if (!cluster_col %in% names(plot_data)) return(data.frame())
  labels <- stats::aggregate(
    cbind(UMAP1, UMAP2) ~ cluster_id,
    data = data.frame(
      cluster_id = as.character(plot_data[[cluster_col]]),
      UMAP1 = suppressWarnings(as.numeric(plot_data$UMAP1)),
      UMAP2 = suppressWarnings(as.numeric(plot_data$UMAP2))
    ),
    FUN = stats::median
  )
  labels[order(labels$cluster_id), , drop = FALSE]
}

add_cluster_outline_layers <- function(plot, clustered_plot_data, label_clusters = TRUE) {
  hulls <- cluster_hull_data(clustered_plot_data)
  if (!nrow(hulls)) return(plot)
  pal <- cluster_palette(hulls$cluster_id)
  for (cluster_id in names(pal)) {
    h <- hulls[hulls$cluster_id == cluster_id, , drop = FALSE]
    plot <- plot +
      ggplot2::geom_path(
        data = h,
        ggplot2::aes(x = UMAP1, y = UMAP2),
        inherit.aes = FALSE,
        color = pal[[cluster_id]],
        linewidth = 0.72,
        linetype = "dashed",
        lineend = "round",
        show.legend = FALSE
      )
  }
  if (isTRUE(label_clusters)) {
    labels <- cluster_label_data(clustered_plot_data)
    for (cluster_id in names(pal)) {
      lab <- labels[labels$cluster_id == cluster_id, , drop = FALSE]
      if (!nrow(lab)) next
      plot <- plot +
        ggplot2::geom_label(
          data = lab,
          ggplot2::aes(x = UMAP1, y = UMAP2, label = cluster_id),
          inherit.aes = FALSE,
          color = pal[[cluster_id]],
          fill = "white",
          linewidth = 0.18,
          size = 2.4,
          fontface = "bold",
          show.legend = FALSE
        )
    }
  }
  plot
}

build_umap_cluster_plot <- function(clustered_plot_data,
                                    initial_size = 0.22,
                                    best_size = 1.2,
                                    shape_by_pred = TRUE) {
  add_cluster_outline_layers(
    build_umap_plot(
      clustered_plot_data,
      initial_size = initial_size,
      best_size = best_size,
      shape_by_pred = shape_by_pred
    ),
    clustered_plot_data
  )
}

save_plot_pair <- function(plot, output_prefix, width = 6.5, height = 6.5) {
  ggplot2::ggsave(
    paste0(output_prefix, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf"
  )
  ggplot2::ggsave(
    paste0(output_prefix, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  message("Wrote PDF: ", paste0(output_prefix, ".pdf"))
  message("Wrote PNG: ", paste0(output_prefix, ".png"))
}

write_umap_outputs <- function(plot_data,
                               figure_prefix,
                               table_prefix,
                               initial_size = 0.22,
                               best_size = 1.2,
                               shape_by_pred = TRUE,
                               cluster_feature_mat = NULL,
                               figures_wclusters_dir = NULL,
                               tables_wclusters_dir = NULL,
                               cluster_seed = 123L,
                               cluster_k_min = 2L,
                               cluster_k_max = 8L,
                               cluster_silhouette_sample_n = 5000L) {
  save_plot_pair(
    build_umap_plot(
      plot_data,
      initial_size = initial_size,
      best_size = best_size,
      shape_by_pred = shape_by_pred
    ),
    figure_prefix
  )
  coord_path <- paste0(table_prefix, "_coordinates.csv")
  utils::write.csv(plot_data, coord_path, quote = FALSE, row.names = FALSE)
  message("Wrote coordinates: ", coord_path)

  if (!is.null(figures_wclusters_dir) || !is.null(tables_wclusters_dir)) {
    write_umap_cluster_outputs(
      plot_data = plot_data,
      feature_mat = cluster_feature_mat,
      figure_prefix = figure_prefix,
      table_prefix = table_prefix,
      initial_size = initial_size,
      best_size = best_size,
      shape_by_pred = shape_by_pred,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )
  }
}

write_umap_cluster_outputs <- function(plot_data,
                                       feature_mat,
                                       figure_prefix,
                                       table_prefix,
                                       initial_size = 0.22,
                                       best_size = 1.2,
                                       shape_by_pred = TRUE,
                                       figures_wclusters_dir,
                                       tables_wclusters_dir,
                                       cluster_seed = 123L,
                                       cluster_k_min = 2L,
                                       cluster_k_max = 8L,
                                       cluster_silhouette_sample_n = 5000L) {
  if (is.null(feature_mat)) {
    stop("cluster_feature_mat must be supplied when Wclusters output is requested.")
  }
  figures_wclusters_dir <- normalizePath(path.expand(figures_wclusters_dir), mustWork = FALSE)
  tables_wclusters_dir <- normalizePath(path.expand(tables_wclusters_dir), mustWork = FALSE)
  subdirs <- cluster_output_subdirs()
  basis <- list(
    umap_coordinates = as.matrix(plot_data[, c("UMAP1", "UMAP2"), drop = FALSE]),
    input_features = feature_mat
  )
  source_labels <- c(
    umap_coordinates = "UMAP1_UMAP2",
    input_features = "input_features"
  )

  for (source_name in names(basis)) {
    source_dir <- subdirs[[source_name]]
    cluster_assignment <- auto_silhouette_kmeans(
      basis_mat = basis[[source_name]],
      plot_data = plot_data,
      cluster_source = source_labels[[source_name]],
      seed = as.integer(cluster_seed) + match(source_name, names(basis)) - 1L,
      k_min = cluster_k_min,
      k_max = cluster_k_max,
      silhouette_sample_n = cluster_silhouette_sample_n
    )
    clustered_plot_data <- plot_data
    clustered_plot_data$cluster_source <- source_labels[[source_name]]
    clustered_plot_data$cluster_id <- cluster_assignment$cluster_id
    clustered_plot_data$cluster_num <- cluster_assignment$cluster_num
    clustered_plot_data$cluster_k <- length(unique(cluster_assignment$cluster_num))
    selected_summary <- cluster_assignment$summary[cluster_assignment$summary$selected, , drop = FALSE]
    clustered_plot_data$cluster_silhouette_avg <- selected_summary$average_silhouette[[1]]
    clustered_plot_data$cluster_silhouette_sample_n <- selected_summary$sample_n[[1]]

    clustered_figure_prefix <- file.path(figures_wclusters_dir, source_dir, basename(figure_prefix))
    clustered_table_prefix <- file.path(tables_wclusters_dir, source_dir, basename(table_prefix))
    dir.create(dirname(clustered_figure_prefix), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(clustered_table_prefix), recursive = TRUE, showWarnings = FALSE)

    save_plot_pair(
      build_umap_cluster_plot(
        clustered_plot_data,
        initial_size = initial_size,
        best_size = best_size,
        shape_by_pred = shape_by_pred
      ),
      clustered_figure_prefix
    )
    coord_path <- paste0(clustered_table_prefix, "_coordinates.csv")
    silhouette_path <- paste0(clustered_table_prefix, "_silhouette.csv")
    utils::write.csv(clustered_plot_data, coord_path, quote = FALSE, row.names = FALSE)
    utils::write.csv(cluster_assignment$summary, silhouette_path, quote = FALSE, row.names = FALSE)
    message("Wrote clustered coordinates: ", coord_path)
    message("Wrote silhouette summary: ", silhouette_path)
  }
}

plot_pca_variance <- function(variance_df, figure_prefix) {
  variance_df$PC_count <- seq_len(nrow(variance_df))
  variance_df$PC_count_label <- factor(variance_df$PC_count, levels = variance_df$PC_count)
  cumulative_plot <- ggplot2::ggplot(
    variance_df,
    ggplot2::aes(x = PC_count_label, y = cumulative_variance_fraction)
  ) +
    ggplot2::geom_col(fill = "#4C78A8", width = 0.72) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(x = "Number of PCs", y = "Cumulative explained variance") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
  save_plot_pair(cumulative_plot, figure_prefix, width = 6.2, height = 4.4)

  variance_df$PC <- factor(variance_df$PC, levels = variance_df$PC)
  individual_plot <- ggplot2::ggplot(
    variance_df,
    ggplot2::aes(x = PC, y = variance_fraction)
  ) +
    ggplot2::geom_col(fill = "#4C78A8", width = 0.72) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x, 1), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(x = "PC", y = "Explained variance") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
  save_plot_pair(individual_plot, paste0(figure_prefix, "_individual"), width = 6.2, height = 4.4)
}

run_pca <- function(feature_mat, pca_n, variance_path, variance_figure_prefix) {
  pca_n <- as.integer(pca_n)
  if (!is.finite(pca_n) || is.na(pca_n) || pca_n < 1L) {
    stop("pca_n must be a positive integer.")
  }
  pca_n <- min(pca_n, ncol(feature_mat))
  message("Running PCA on standardized UMAP features; retaining ", pca_n, " PCs.")
  pca <- stats::prcomp(feature_mat, center = FALSE, scale. = FALSE)
  variance <- pca$sdev^2
  variance_df <- data.frame(
    PC = paste0("PC", seq_along(variance)),
    variance = variance,
    variance_fraction = variance / sum(variance),
    cumulative_variance_fraction = cumsum(variance) / sum(variance)
  )
  utils::write.csv(variance_df, variance_path, quote = FALSE, row.names = FALSE)
  message("Wrote PCA variance: ", variance_path)
  plot_pca_variance(variance_df, variance_figure_prefix)
  pca$x[, seq_len(pca_n), drop = FALSE]
}

sample_initial_rows <- function(initial_df, sample_n, seed, by_seed = FALSE) {
  sample_n <- as.integer(sample_n)
  if (!is.finite(sample_n) || is.na(sample_n) || sample_n < 1L) {
    stop("sample_initial_n must be a positive integer.")
  }
  if (isTRUE(by_seed)) {
    if (!"seed" %in% names(initial_df)) stop("Initial population CSV must contain a seed column.")
    split_rows <- split(seq_len(nrow(initial_df)), as.integer(initial_df$seed))
    seed_ids <- sort(as.integer(names(split_rows)))
    set.seed(as.integer(seed))
    sampled <- vapply(
      seed_ids,
      function(seed_id) sample(split_rows[[as.character(seed_id)]], size = 1L),
      integer(1)
    )
    message("Sampled one initial point per seed across ", length(sampled), " seeds.")
    return(as.integer(sampled))
  }
  if (sample_n >= nrow(initial_df)) {
    message("Requested sampled initial UMAP with ", sample_n, " initial points; using all ", nrow(initial_df), " points.")
    return(seq_len(nrow(initial_df)))
  }
  set.seed(as.integer(seed))
  sort(sample.int(nrow(initial_df), size = sample_n, replace = FALSE))
}

default_invivo_data_dir <- function() {
  oxygen_dir <- default_oxygen_dir()
  candidates <- c(
    file.path(oxygen_dir, "..", "data", "InVivoData_Gemcitabine"),
    file.path(oxygen_dir, "data", "InVivoData_Gemcitabine"),
    file.path(dirname(oxygen_dir), "data", "InVivoData_Gemcitabine")
  )
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "dt_Gem_VT_20260209_v5.xlsx")) &&
        file.exists(file.path(candidate, "all_ploidy.csv"))) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  normalizePath(candidates[[1L]], mustWork = FALSE)
}

invivo_data_paths <- function(data_dir = default_invivo_data_dir()) {
  data_dir <- normalizePath(path.expand(data_dir), mustWork = FALSE)
  dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
  ploidy_path <- file.path(data_dir, "all_ploidy.csv")
  if (!file.exists(dt_path)) stop("Missing in vivo tumor-burden file: ", dt_path)
  if (!file.exists(ploidy_path)) stop("Missing in vivo ploidy file: ", ploidy_path)
  list(dt_path = dt_path, ploidy_path = ploidy_path)
}

invivo_simulation_env <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) return(cache)
    workflow_root <- normalizePath(file.path(default_oxygen_dir(), "code", "O2_supply_demand_MAP"), mustWork = FALSE)
    shared_path <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
    common_path <- file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R")
    model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
    viz_path <- file.path(workflow_root, "vis", "viz_invivo_model_O2_supply_demand_MAP_results.R")
    if (!file.exists(shared_path)) stop("Shared utility script not found: ", shared_path)
    if (!file.exists(common_path)) stop("Common semantics script not found: ", common_path)
    if (!file.exists(model_path)) stop("Model script not found: ", model_path)
    if (!file.exists(viz_path)) stop("In vivo visualization script not found: ", viz_path)

    env <- new.env(parent = globalenv())
    env$commandArgs <- function(trailingOnly = FALSE) character(0)
    sys.source(shared_path, envir = globalenv(), chdir = TRUE)
    sys.source(common_path, envir = globalenv(), chdir = TRUE)
    sys.source(shared_path, envir = env, chdir = TRUE)
    sys.source(common_path, envir = env, chdir = TRUE)
    Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
    sys.source(model_path, envir = env, chdir = TRUE)
    sys.source(viz_path, envir = env, chdir = TRUE)

    required <- c(
      "normalize_cfg_for_viz",
      "read_run_params",
      "prepare_data",
      "simulate_one_full_horizon",
      ".lambda_eff_of_O2",
      ".mu_eff_of_O2",
      ".pmisseg_of_O2",
      ".p_wgd_of_O2",
      ".pr_delta_vec"
    )
    missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
    if (length(missing)) stop("In vivo simulation environment is missing: ", paste(missing, collapse = ", "))
    cache <<- env
    cache
  }
})

compute_misseg_death_rate <- function(lambda,
                                      O2,
                                      N,
                                      run_params,
                                      cfg = NULL,
                                      model_env = invivo_simulation_env()) {
  lambda <- suppressWarnings(as.numeric(lambda))
  O2 <- suppressWarnings(as.numeric(O2))
  N <- suppressWarnings(as.numeric(N))
  if (!(length(lambda) == length(O2) && length(O2) == length(N))) {
    stop("lambda, O2, and N must have equal lengths.")
  }

  p_mis <- get(".pmisseg_of_O2", envir = model_env, inherits = TRUE)(
    O2 = O2,
    run_params = run_params,
    N = N
  )
  p_wgd <- get(".p_wgd_of_O2", envir = model_env, inherits = TRUE)(
    O2 = O2,
    run_params = run_params
  )
  pr_delta <- get(".pr_delta_vec", envir = model_env, inherits = TRUE)

  buffer_smax <- as.numeric(run_params$buffer_smax %||% 1.0)
  if (!is.finite(buffer_smax)) buffer_smax <- 1.0
  buffer_beta <- as.numeric(run_params$buffer_beta %||% 0.0)
  if (!is.finite(buffer_beta)) buffer_beta <- 0.0
  buffer_n_exp <- as.numeric(run_params$buffer_n_exp %||% 1.0)
  if (!is.finite(buffer_n_exp)) buffer_n_exp <- 1.0
  n_unit <- as.integer((cfg$N_UNIT %||% cfg$N_unit %||% run_params$N_UNIT %||% 22L)[[1L]])
  if (!is.finite(n_unit) || n_unit <= 0L) n_unit <- 22L

  cache <- new.env(parent = emptyenv())
  mass_dropped <- vapply(seq_along(N), function(i) {
    if (!is.finite(N[[i]]) || !is.finite(p_mis[[i]])) return(NA_real_)
    n_i <- as.integer(round(N[[i]]))
    p_i <- pmin(pmax(as.numeric(p_mis[[i]]), 0), 1)
    key <- paste(n_i, signif(p_i, 12), sep = "|")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    delta <- pr_delta(
      N = n_i,
      p = p_i,
      eps_tail = 1e-8,
      N_unit = n_unit,
      buffer_smax = buffer_smax,
      buffer_beta = buffer_beta,
      buffer_n_exp = buffer_n_exp
    )
    value <- as.numeric(attr(delta, "mass_dropped"))
    if (!is.finite(value)) value <- NA_real_
    assign(key, value, envir = cache)
    value
  }, numeric(1))

  dead_daughters_per_division <- pmin(pmax(2 * mass_dropped, 0), 2)
  p_wgd <- pmin(pmax(as.numeric(p_wgd), 0), 1)
  lambda * (1 - p_wgd) * dead_daughters_per_division
}

time_weighted_mean <- function(day, value) {
  day <- suppressWarnings(as.numeric(day))
  value <- suppressWarnings(as.numeric(value))
  keep <- is.finite(day) & is.finite(value)
  day <- day[keep]
  value <- value[keep]
  if (!length(day)) return(NA_real_)

  ord <- order(day)
  day <- day[ord]
  value <- value[ord]
  day_unique <- sort(unique(day))
  if (length(day_unique) < length(day)) {
    value <- vapply(day_unique, function(d) mean(value[day == d], na.rm = TRUE), numeric(1))
    day <- day_unique
  }
  if (length(day) == 1L || max(day) <= min(day)) return(mean(value, na.rm = TRUE))

  widths <- diff(day)
  integral <- sum(widths * (head(value, -1L) + tail(value, -1L)) / 2, na.rm = TRUE)
  integral / (max(day) - min(day))
}

compute_invivo_rate_summary_for_sim <- function(sim,
                                                run_params,
                                                cfg,
                                                model_env,
                                                horizon_day = 100) {
  burden <- sim$burden
  ploidy <- sim$ploidy
  if (is.null(burden) || is.null(ploidy) || !nrow(burden) || !nrow(ploidy)) {
    return(NULL)
  }
  burden <- burden[is.finite(burden$day) & burden$day <= horizon_day + 1e-9, , drop = FALSE]
  ploidy <- ploidy[is.finite(ploidy$day) & ploidy$day <= horizon_day + 1e-9, , drop = FALSE]
  if (!nrow(burden) || !nrow(ploidy)) return(NULL)
  mean_o2 <- time_weighted_mean(burden$day, burden$pred_o2_pct)

  o2_by_day <- stats::setNames(as.numeric(burden$pred_o2_pct), as.character(as.numeric(burden$day)))
  ploidy$o2_pct <- as.numeric(o2_by_day[as.character(as.numeric(ploidy$day))])
  keep <- is.finite(ploidy$o2_pct) & is.finite(ploidy$N) & is.finite(ploidy$fraction)
  ploidy <- ploidy[keep, , drop = FALSE]
  if (!nrow(ploidy)) return(NULL)

  lambda <- get(".lambda_eff_of_O2", envir = model_env, inherits = TRUE)(
    O2 = ploidy$o2_pct,
    run_params = run_params,
    N = ploidy$N,
    O2_growth = isTRUE(cfg$O2_growth %||% TRUE)
  )
  mu <- get(".mu_eff_of_O2", envir = model_env, inherits = TRUE)(
    O2 = ploidy$o2_pct,
    run_params = run_params,
    N = ploidy$N
  )
  misseg_death <- compute_misseg_death_rate(
    lambda = lambda,
    O2 = ploidy$o2_pct,
    N = ploidy$N,
    run_params = run_params,
    cfg = cfg,
    model_env = model_env
  )

  rate_df <- data.frame(
    day = as.numeric(ploidy$day),
    fraction = pmax(as.numeric(ploidy$fraction), 0),
    lambda = as.numeric(lambda),
    mu = as.numeric(mu),
    misseg_death = as.numeric(misseg_death),
    stringsAsFactors = FALSE
  )
  rate_df <- rate_df[is.finite(rate_df$day) & is.finite(rate_df$fraction) &
                       is.finite(rate_df$lambda) & is.finite(rate_df$mu) &
                       is.finite(rate_df$misseg_death), , drop = FALSE]
  if (!nrow(rate_df)) return(NULL)

  by_day <- split(rate_df, rate_df$day)
  day_rates <- do.call(rbind, lapply(by_day, function(df) {
    w <- pmax(df$fraction, 0)
    w_sum <- sum(w, na.rm = TRUE)
    if (!is.finite(w_sum) || w_sum <= 0) return(NULL)
    lambda_mean <- sum(w * df$lambda, na.rm = TRUE) / w_sum
    mu_mean <- sum(w * df$mu, na.rm = TRUE) / w_sum
    misseg_death_mean <- sum(w * df$misseg_death, na.rm = TRUE) / w_sum
    data.frame(
      day = as.numeric(df$day[[1L]]),
      net_growth = lambda_mean - mu_mean,
      turnover_rate = lambda_mean + mu_mean,
      misseg_death_rate = misseg_death_mean,
      net_growth_with_misseg_death = lambda_mean - mu_mean - misseg_death_mean,
      turnover_rate_with_misseg_death = lambda_mean + mu_mean + misseg_death_mean,
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(day_rates) || !nrow(day_rates)) return(NULL)

  data.frame(
    mean_net_growth_0_100d = time_weighted_mean(day_rates$day, day_rates$net_growth),
    mean_turnover_rate_0_100d = time_weighted_mean(day_rates$day, day_rates$turnover_rate),
    mean_O2_0_100d = mean_o2,
    mean_misseg_death_rate_0_100d = time_weighted_mean(day_rates$day, day_rates$misseg_death_rate),
    mean_net_growth_with_misseg_death_0_100d = time_weighted_mean(day_rates$day, day_rates$net_growth_with_misseg_death),
    mean_turnover_rate_with_misseg_death_0_100d = time_weighted_mean(day_rates$day, day_rates$turnover_rate_with_misseg_death),
    stringsAsFactors = FALSE
  )
}

select_invivo_representative_scenarios <- function(scenarios) {
  cohorts <- vapply(scenarios, function(sc) as.character(sc$cohort %||% ""), character(1))
  reps <- lapply(c("2N", "4N"), function(cohort_id) {
    idx <- which(cohorts == cohort_id)
    if (!length(idx)) stop("No ", cohort_id, " scenario found for growth/turnover calculation.")
    scenarios[[idx[[1L]]]]
  })
  names(reps) <- c("2N", "4N")
  reps
}

compute_invivo_growth_turnover_one_seed <- function(seed_dir,
                                                    data_paths,
                                                    model_env = invivo_simulation_env(),
                                                    horizon_day = 100,
                                                    report_dt = 1.0) {
  seed <- seed_from_dir(seed_dir)
  cfg_path <- file.path(seed_dir, "fit_config.rds")
  if (!file.exists(cfg_path)) stop("Missing fit_config.rds: ", cfg_path)
  cfg <- readRDS(cfg_path)
  cfg <- get("normalize_cfg_for_viz", envir = model_env, inherits = TRUE)(cfg)
  run_params <- get("read_run_params", envir = model_env, inherits = TRUE)(seed_dir, cfg = cfg)
  scenarios <- get("prepare_data", envir = model_env, inherits = TRUE)(
    data_paths$dt_path,
    data_paths$ploidy_path,
    cfg
  )
  scenarios <- select_invivo_representative_scenarios(scenarios)

  sim_fun <- get("simulate_one_full_horizon", envir = model_env, inherits = TRUE)
  scenario_rows <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    sim <- sim_fun(
      run_params = run_params,
      scenario = sc,
      cfg = cfg,
      horizon_day = horizon_day,
      report_dt = report_dt
    )
    rates <- compute_invivo_rate_summary_for_sim(
      sim = sim,
      run_params = run_params,
      cfg = cfg,
      model_env = model_env,
      horizon_day = horizon_day
    )
    if (is.null(rates)) next
    scenario_rows[[i]] <- data.frame(
      seed = seed,
      cohort = as.character(sc$cohort),
      harvest = as.character(sc$harvest),
      dose = as.numeric(sc$dose),
      rates,
      stringsAsFactors = FALSE
    )
  }
  scenario_df <- do.call(rbind, scenario_rows)
  if (is.null(scenario_df) || !nrow(scenario_df)) {
    stop("No rate summaries were generated for seed", seed)
  }

  cohort_levels <- c("2N", "4N")
  cohort_summary <- do.call(rbind, lapply(cohort_levels, function(cohort_id) {
    df <- scenario_df[scenario_df$cohort == cohort_id, , drop = FALSE]
    if (!nrow(df)) return(NULL)
    data.frame(
      cohort = cohort_id,
      mean_net_growth_0_100d = mean(df$mean_net_growth_0_100d, na.rm = TRUE),
      mean_turnover_rate_0_100d = mean(df$mean_turnover_rate_0_100d, na.rm = TRUE),
      mean_O2_0_100d = mean(df$mean_O2_0_100d, na.rm = TRUE),
      mean_misseg_death_rate_0_100d = mean(df$mean_misseg_death_rate_0_100d, na.rm = TRUE),
      mean_net_growth_with_misseg_death_0_100d = mean(df$mean_net_growth_with_misseg_death_0_100d, na.rm = TRUE),
      mean_turnover_rate_with_misseg_death_0_100d = mean(df$mean_turnover_rate_with_misseg_death_0_100d, na.rm = TRUE),
      n_scenarios = nrow(df),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(cohort_summary) || !nrow(cohort_summary)) {
    stop("No 2N/4N cohort summaries were generated for seed", seed)
  }

  value_for <- function(cohort_id, col) {
    hit <- cohort_summary[cohort_summary$cohort == cohort_id, col, drop = TRUE]
    if (length(hit)) as.numeric(hit[[1L]]) else NA_real_
  }
  n_for <- function(cohort_id) {
    hit <- cohort_summary[cohort_summary$cohort == cohort_id, "n_scenarios", drop = TRUE]
    if (length(hit)) as.integer(hit[[1L]]) else 0L
  }

  net_2n <- value_for("2N", "mean_net_growth_0_100d")
  net_4n <- value_for("4N", "mean_net_growth_0_100d")
  turnover_2n <- value_for("2N", "mean_turnover_rate_0_100d")
  turnover_4n <- value_for("4N", "mean_turnover_rate_0_100d")
  o2_2n <- value_for("2N", "mean_O2_0_100d")
  o2_4n <- value_for("4N", "mean_O2_0_100d")
  misseg_death_2n <- value_for("2N", "mean_misseg_death_rate_0_100d")
  misseg_death_4n <- value_for("4N", "mean_misseg_death_rate_0_100d")
  net_with_misseg_death_2n <- value_for("2N", "mean_net_growth_with_misseg_death_0_100d")
  net_with_misseg_death_4n <- value_for("4N", "mean_net_growth_with_misseg_death_0_100d")
  turnover_with_misseg_death_2n <- value_for("2N", "mean_turnover_rate_with_misseg_death_0_100d")
  turnover_with_misseg_death_4n <- value_for("4N", "mean_turnover_rate_with_misseg_death_0_100d")

  data.frame(
    seed = seed,
    mean_net_growth_0_100d_2N = net_2n,
    mean_net_growth_0_100d_4N = net_4n,
    mean_net_growth_0_100d = mean(c(net_2n, net_4n), na.rm = TRUE),
    mean_turnover_rate_0_100d_2N = turnover_2n,
    mean_turnover_rate_0_100d_4N = turnover_4n,
    mean_turnover_rate_0_100d = mean(c(turnover_2n, turnover_4n), na.rm = TRUE),
    mean_O2_0_100d_2N = o2_2n,
    mean_O2_0_100d_4N = o2_4n,
    mean_O2_0_100d = mean(c(o2_2n, o2_4n), na.rm = TRUE),
    mean_misseg_death_rate_0_100d_2N = misseg_death_2n,
    mean_misseg_death_rate_0_100d_4N = misseg_death_4n,
    mean_misseg_death_rate_0_100d = mean(c(misseg_death_2n, misseg_death_4n), na.rm = TRUE),
    mean_net_growth_with_misseg_death_0_100d_2N = net_with_misseg_death_2n,
    mean_net_growth_with_misseg_death_0_100d_4N = net_with_misseg_death_4n,
    mean_net_growth_with_misseg_death_0_100d = mean(c(net_with_misseg_death_2n, net_with_misseg_death_4n), na.rm = TRUE),
    mean_turnover_rate_with_misseg_death_0_100d_2N = turnover_with_misseg_death_2n,
    mean_turnover_rate_with_misseg_death_0_100d_4N = turnover_with_misseg_death_4n,
    mean_turnover_rate_with_misseg_death_0_100d = mean(c(turnover_with_misseg_death_2n, turnover_with_misseg_death_4n), na.rm = TRUE),
    n_rate_scenarios_2N = n_for("2N"),
    n_rate_scenarios_4N = n_for("4N"),
    n_rate_scenarios = nrow(scenario_df),
    rate_horizon_day = as.numeric(horizon_day),
    rate_report_dt = as.numeric(report_dt),
    stringsAsFactors = FALSE
  )
}

paper_generate_invivo_growth_turnover_table <- function(input_dir = default_dataset_input_dir("invivo"),
                                                        best_csv = file.path(paper_tables_dir("invivo"), "invivo_best_params_by_seed.csv"),
                                                        output_csv = file.path(paper_tables_dir("invivo"), "invivo_best_params_growth_turnover_100d.csv"),
                                                        data_dir = default_invivo_data_dir(),
                                                        max_seeds = NA_integer_,
                                                        horizon_day = 100,
                                                        report_dt = 1.0) {
  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  best_csv <- normalizePath(path.expand(best_csv), mustWork = FALSE)
  output_csv <- normalizePath(path.expand(output_csv), mustWork = FALSE)
  max_seeds <- as_int(max_seeds, NA_integer_)
  horizon_day <- as_num(horizon_day, 100)
  report_dt <- as_num(report_dt, 1.0)
  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  if (!file.exists(best_csv)) stop("Missing best-parameter CSV: ", best_csv)
  if (!is.finite(horizon_day) || horizon_day <= 0) stop("horizon_day must be > 0.")
  if (!is.finite(report_dt) || report_dt <= 0) stop("report_dt must be > 0.")

  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  best_df <- read_csv_plain(best_csv)
  if (!"seed" %in% names(best_df)) stop("Best-parameter CSV must contain a seed column.")
  best_df$seed <- as.integer(best_df$seed)

  seed_dirs <- list_seed_dirs(input_dir)
  if (!length(seed_dirs)) stop("No seed directories found under: ", input_dir)
  if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
    seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
  }
  seed_ids <- vapply(seed_dirs, seed_from_dir, integer(1))
  missing_best <- setdiff(seed_ids, best_df$seed)
  if (length(missing_best)) {
    stop("Best-parameter CSV is missing seeds: ", paste(head(missing_best, 20), collapse = ", "))
  }

  data_paths <- invivo_data_paths(data_dir)
  model_env <- invivo_simulation_env()
  message("Computing in vivo 0-", horizon_day, " day growth/turnover metrics for ", length(seed_dirs), " seeds.")
  rows <- vector("list", length(seed_dirs))
  for (i in seq_along(seed_dirs)) {
    rows[[i]] <- compute_invivo_growth_turnover_one_seed(
      seed_dir = seed_dirs[[i]],
      data_paths = data_paths,
      model_env = model_env,
      horizon_day = horizon_day,
      report_dt = report_dt
    )
    if (i %% 10L == 0L || i == length(seed_dirs)) {
      message("Computed growth/turnover metrics for ", i, "/", length(seed_dirs), " seeds.")
    }
  }
  metrics_df <- do.call(rbind, rows)
  best_subset <- best_df[match(metrics_df$seed, best_df$seed), , drop = FALSE]
  metric_cols <- setdiff(names(metrics_df), "seed")
  out <- cbind(best_subset, metrics_df[, metric_cols, drop = FALSE])
  out <- append_invivo_pred1000_ploidy_ratios(out, input_dir = input_dir)
  required_rate_cols <- c(
    "mean_net_growth_0_100d",
    "mean_turnover_rate_0_100d",
    "mean_O2_0_100d",
    "mean_misseg_death_rate_0_100d",
    "mean_net_growth_with_misseg_death_0_100d",
    "mean_turnover_rate_with_misseg_death_0_100d"
  )
  if (any(vapply(required_rate_cols, function(col) any(!is.finite(out[[col]])), logical(1)))) {
    stop("Computed growth/turnover table contains non-finite aggregate rate metrics.")
  }
  utils::write.csv(out, output_csv, quote = FALSE, row.names = FALSE)
  message("Wrote in vivo growth/turnover table: ", output_csv)
  invisible(output_csv)
}

append_invivo_pred1000_ploidy_ratios <- function(df,
                                                 input_dir = default_dataset_input_dir("invivo"),
                                                 target_day = 1000) {
  if (!"seed" %in% names(df)) stop("Input table must contain a seed column.")
  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  seed_ids <- as.integer(df$seed)
  rows <- lapply(seed_ids, function(seed) {
    seed_dir <- file.path(input_dir, paste0("seed", seed))
    read_pred1000_ploidy_ratios(seed_dir, seed, target_day = target_day)
  })
  ratio_df <- do.call(rbind, rows)
  for (col in names(ratio_df)) {
    df[[col]] <- ratio_df[[col]]
  }
  ratio_cols <- names(ratio_df)
  if (any(vapply(ratio_cols, function(col) any(!is.finite(df[[col]])), logical(1)))) {
    stop("Computed pred1000 ploidy ratio columns contain non-finite values.")
  }
  df
}

paper_update_invivo_growth_turnover_ploidy_ratios <- function(input_dir = default_dataset_input_dir("invivo"),
                                                              growth_turnover_csv = file.path(paper_tables_dir("invivo"), "invivo_best_params_growth_turnover_100d.csv"),
                                                              target_day = 1000) {
  input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
  growth_turnover_csv <- normalizePath(path.expand(growth_turnover_csv), mustWork = FALSE)
  if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
  if (!file.exists(growth_turnover_csv)) stop("Missing growth/turnover CSV: ", growth_turnover_csv)
  df <- read_csv_plain(growth_turnover_csv)
  out <- append_invivo_pred1000_ploidy_ratios(df, input_dir = input_dir, target_day = target_day)
  utils::write.csv(out, growth_turnover_csv, quote = FALSE, row.names = FALSE)
  message("Updated in vivo growth/turnover table with pred1000 ploidy ratios: ", growth_turnover_csv)
  invisible(growth_turnover_csv)
}

transform_extra_numeric_features <- function(df, feature_cols) {
  feature_cols <- as.character(feature_cols)
  missing <- setdiff(feature_cols, names(df))
  if (length(missing)) stop("Input table is missing extra feature columns: ", paste(missing, collapse = ", "))
  out <- as.data.frame(lapply(feature_cols, function(col) {
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (any(!is.finite(vals))) stop("Extra UMAP feature contains non-finite values: ", col)
    vals
  }), check.names = FALSE)
  names(out) <- feature_cols
  out
}

paper_generate_invivo_growth_turnover_umap_figures <- function(tables_dir = paper_tables_dir("invivo"),
                                                               figures_dir = paper_figures_dir("invivo"),
                                                               support_tables_dir = tables_dir,
                                                               figures_wclusters_dir = NULL,
                                                               tables_wclusters_dir = NULL,
                                                               growth_turnover_csv = file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
                                                               objective_seed_dir = default_dataset_input_dir("invivo"),
                                                               output_prefix = "invivo_best_params_growth_turnover_100d_umap",
                                                               run_clustered_umaps = FALSE,
                                                               shape_by_pred = TRUE,
                                                               umap_seed = 123L,
                                                               n_neighbors = 80L,
                                                               min_dist = 0.1,
                                                               n_threads = max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L)),
                                                               pca_n = 10L,
                                                               best_size = 1.6,
                                                               cluster_seed = 123L,
                                                               cluster_k_min = 2L,
                                                               cluster_k_max = 8L,
                                                               cluster_silhouette_sample_n = 5000L,
                                                               rate_version = c(
                                                                 "base",
                                                                 "with_misseg_death",
                                                                 "base_with_o2",
                                                                 "with_misseg_death_with_o2",
                                                                 "base_with_o2_ploidy_ratio",
                                                                 "with_misseg_death_with_o2_ploidy_ratio"
                                                               )) {
  for (pkg in c("ggplot2", "uwot")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required R package is not installed: ", pkg)
    }
  }
  rate_version <- match.arg(rate_version)
  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  figures_dir <- normalizePath(path.expand(figures_dir), mustWork = FALSE)
  support_tables_dir <- normalizePath(path.expand(support_tables_dir), mustWork = FALSE)
  growth_turnover_csv <- normalizePath(path.expand(growth_turnover_csv), mustWork = FALSE)
  objective_seed_dir <- normalizePath(path.expand(objective_seed_dir), mustWork = FALSE)
  if (isTRUE(run_clustered_umaps)) {
    figures_wclusters_dir <- normalizePath(path.expand(figures_wclusters_dir %||% file.path(dirname(figures_dir), "FiguresWclusters")), mustWork = FALSE)
    tables_wclusters_dir <- normalizePath(path.expand(tables_wclusters_dir %||% file.path(dirname(support_tables_dir), "TablesWclusters")), mustWork = FALSE)
  } else {
    figures_wclusters_dir <- NULL
    tables_wclusters_dir <- NULL
  }
  if (!file.exists(growth_turnover_csv)) stop("Missing growth/turnover CSV: ", growth_turnover_csv)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(support_tables_dir, recursive = TRUE, showWarnings = FALSE)

  params <- umap_parameter_set("invivo")
  log10_params <- umap_log10_parameter_set("invivo")
  extra_features <- if (identical(rate_version, "with_misseg_death")) {
    c("mean_net_growth_with_misseg_death_0_100d", "mean_turnover_rate_with_misseg_death_0_100d")
  } else if (identical(rate_version, "with_misseg_death_with_o2")) {
    c("mean_net_growth_with_misseg_death_0_100d", "mean_turnover_rate_with_misseg_death_0_100d", "mean_O2_0_100d")
  } else if (identical(rate_version, "with_misseg_death_with_o2_ploidy_ratio")) {
    c(
      "mean_net_growth_with_misseg_death_0_100d",
      "mean_turnover_rate_with_misseg_death_0_100d",
      "mean_O2_0_100d",
      "pred1000_ploidy_ratio_2N",
      "pred1000_ploidy_ratio_4N"
    )
  } else if (identical(rate_version, "base_with_o2_ploidy_ratio")) {
    c(
      "mean_net_growth_0_100d",
      "mean_turnover_rate_0_100d",
      "mean_O2_0_100d",
      "pred1000_ploidy_ratio_2N",
      "pred1000_ploidy_ratio_4N"
    )
  } else if (identical(rate_version, "base_with_o2")) {
    c("mean_net_growth_0_100d", "mean_turnover_rate_0_100d", "mean_O2_0_100d")
  } else {
    c("mean_net_growth_0_100d", "mean_turnover_rate_0_100d")
  }
  reduction_prefix <- if (identical(rate_version, "with_misseg_death")) {
    "best_only_growth_turnover_with_misseg_death"
  } else if (identical(rate_version, "with_misseg_death_with_o2")) {
    "best_only_growth_turnover_with_misseg_death_o2"
  } else if (identical(rate_version, "with_misseg_death_with_o2_ploidy_ratio")) {
    "best_only_growth_turnover_with_misseg_death_o2_ploidy_ratio"
  } else if (identical(rate_version, "base_with_o2_ploidy_ratio")) {
    "best_only_growth_turnover_o2_ploidy_ratio"
  } else if (identical(rate_version, "base_with_o2")) {
    "best_only_growth_turnover_o2"
  } else {
    "best_only_growth_turnover"
  }
  run_label_suffix <- if (identical(rate_version, "with_misseg_death")) {
    " with misseg death"
  } else if (identical(rate_version, "with_misseg_death_with_o2")) {
    " with misseg death and O2"
  } else if (identical(rate_version, "with_misseg_death_with_o2_ploidy_ratio")) {
    " with misseg death, O2, and ploidy ratio"
  } else if (identical(rate_version, "base_with_o2_ploidy_ratio")) {
    " with O2 and ploidy ratio"
  } else if (identical(rate_version, "base_with_o2")) {
    " with O2"
  } else {
    ""
  }
  message("Reading in vivo growth/turnover UMAP input: ", growth_turnover_csv)
  best_df <- attach_objective(read_csv_plain(growth_turnover_csv), objective_seed_dir = objective_seed_dir)
  if (!"seed" %in% names(best_df)) stop("Growth/turnover CSV must contain a seed column.")
  if (isTRUE(shape_by_pred) && !"pred1000_both_gt44" %in% names(best_df)) {
    stop("Growth/turnover CSV must contain pred1000_both_gt44 when shape_by_pred=TRUE.")
  }
  if (isTRUE(shape_by_pred)) {
    best_df$pred1000_both_gt44 <- coerce_logical_column(best_df$pred1000_both_gt44, "pred1000_both_gt44")
  }

  best_features <- cbind(
    transform_umap_features(best_df, params, log10_params),
    transform_extra_numeric_features(best_df, extra_features)
  )
  best_feature_mat <- standardize_features(best_features)
  empty_initial_df <- best_df[0L, , drop = FALSE]

  figure_prefix <- file.path(figures_dir, output_prefix)
  table_prefix <- file.path(support_tables_dir, output_prefix)
  pca_figure_prefix <- default_pca_output_prefix(figure_prefix)
  pca_table_prefix <- default_pca_output_prefix(table_prefix)

  emb <- run_umap_embedding(best_feature_mat, paste0("best-only growth/turnover direct", run_label_suffix), umap_seed, n_neighbors, min_dist, n_threads)
  plot_data <- build_plot_data(
    emb,
    empty_initial_df,
    best_df,
    reduction_label = paste0(reduction_prefix, "_umap"),
    shape_by_pred = shape_by_pred
  )
  write_umap_outputs(
    plot_data,
    figure_prefix,
    table_prefix,
    best_size = best_size,
    shape_by_pred = shape_by_pred,
    cluster_feature_mat = best_feature_mat,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n
  )

  pca_features <- run_pca(
    best_feature_mat,
    pca_n = pca_n,
    variance_path = paste0(pca_table_prefix, "_pca_variance.csv"),
    variance_figure_prefix = paste0(pca_figure_prefix, "_pca_variance_bar")
  )
  pca_emb <- run_umap_embedding(pca_features, paste0("best-only growth/turnover PCA", run_label_suffix), umap_seed, n_neighbors, min_dist, n_threads)
  pca_plot_data <- build_plot_data(
    pca_emb,
    empty_initial_df,
    best_df,
    reduction_label = paste0(reduction_prefix, "_pca", ncol(pca_features), "_umap"),
    shape_by_pred = shape_by_pred
  )
  write_umap_outputs(
    pca_plot_data,
    pca_figure_prefix,
    pca_table_prefix,
    best_size = best_size,
    shape_by_pred = shape_by_pred,
    cluster_feature_mat = pca_features,
    figures_wclusters_dir = figures_wclusters_dir,
    tables_wclusters_dir = tables_wclusters_dir,
    cluster_seed = cluster_seed,
    cluster_k_min = cluster_k_min,
    cluster_k_max = cluster_k_max,
    cluster_silhouette_sample_n = cluster_silhouette_sample_n
  )
  message("Growth/turnover UMAP features: ", paste(c(params, extra_features), collapse = ", "))
  invisible(TRUE)
}

paper_generate_invivo_growth_turnover_boxplot <- function(tables_dir = paper_tables_dir("invivo"),
                                                          figures_dir = paper_figures_dir("invivo"),
                                                          growth_turnover_csv = file.path(tables_dir, "invivo_best_params_growth_turnover_100d.csv"),
                                                          best_csv = file.path(tables_dir, "invivo_best_params_by_seed.csv"),
                                                          output_prefix = "invivo_best_params_growth_turnover_100d_paired_boxplot",
                                                          jitter_seed = 123L,
                                                          jitter_width = 0.22,
                                                          rate_version = c("base", "with_misseg_death")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Required R package is not installed: ggplot2")
  }
  rate_version <- match.arg(rate_version)
  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  figures_dir <- normalizePath(path.expand(figures_dir), mustWork = FALSE)
  growth_turnover_csv <- normalizePath(path.expand(growth_turnover_csv), mustWork = FALSE)
  best_csv <- normalizePath(path.expand(best_csv), mustWork = FALSE)
  if (!file.exists(growth_turnover_csv)) stop("Missing growth/turnover CSV: ", growth_turnover_csv)
  if (!file.exists(best_csv)) stop("Missing best-params CSV: ", best_csv)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  value_cols <- if (identical(rate_version, "with_misseg_death")) {
    c("mean_net_growth_with_misseg_death_0_100d", "mean_turnover_rate_with_misseg_death_0_100d")
  } else {
    c("mean_net_growth_0_100d", "mean_turnover_rate_0_100d")
  }
  required_cols <- c("seed", value_cols)
  df <- read_csv_plain(growth_turnover_csv)
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) stop("Growth/turnover CSV is missing columns: ", paste(missing, collapse = ", "))

  best_df <- read_csv_plain(best_csv)
  best_required_cols <- c("seed", "pred1000_both_gt44")
  best_missing <- setdiff(best_required_cols, names(best_df))
  if (length(best_missing)) stop("Best-params CSV is missing columns: ", paste(best_missing, collapse = ", "))
  if (anyDuplicated(best_df$seed)) stop("Best-params CSV contains duplicated seed values.")
  best_df$pred1000_both_gt44 <- coerce_logical_column(best_df$pred1000_both_gt44, "pred1000_both_gt44")

  seed <- df$seed
  net_growth <- suppressWarnings(as.numeric(df[[value_cols[[1L]]]]))
  turnover <- suppressWarnings(as.numeric(df[[value_cols[[2L]]]]))
  if (any(!is.finite(net_growth))) stop("Non-finite values found in ", value_cols[[1L]], ".")
  if (any(!is.finite(turnover))) stop("Non-finite values found in ", value_cols[[2L]], ".")
  best_idx <- match(seed, best_df$seed)
  if (anyNA(best_idx)) stop("Some growth/turnover seeds are missing from best-params CSV.")
  pred1000_both_gt44 <- best_df$pred1000_both_gt44[best_idx]

  plot_data <- data.frame(
    seed = rep(seed, times = 2L),
    pred1000_both_gt44 = rep(pred1000_both_gt44, times = 2L),
    metric = factor(
      rep(c("Net growth", "Turnover rate"), each = length(seed)),
      levels = c("Net growth", "Turnover rate")
    ),
    value = c(net_growth, turnover),
    stringsAsFactors = FALSE
  )
  set.seed(as.integer(jitter_seed))
  seed_offsets <- stats::runif(length(seed), min = -jitter_width, max = jitter_width)
  plot_data$metric_index <- as.numeric(plot_data$metric)
  plot_data$x_pos <- plot_data$metric_index + rep(seed_offsets, times = 2L)
  plot_data$pred1000_both_gt44_label <- factor(
    as.character(plot_data$pred1000_both_gt44),
    levels = c("FALSE", "TRUE")
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = value)) +
    ggplot2::geom_line(
      ggplot2::aes(x = x_pos, group = seed),
      color = "grey48",
      alpha = 0.16,
      linewidth = 0.22
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(x = metric_index, group = metric),
      width = 0.42,
      outlier.shape = NA,
      fill = "grey78",
      alpha = 0.72,
      linewidth = 0.38,
      color = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = x_pos, fill = pred1000_both_gt44_label, shape = pred1000_both_gt44_label),
      size = 1.05,
      stroke = 0.18,
      color = "black",
      alpha = 0.62
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c("Net growth", "Turnover rate"),
      limits = c(0.58, 2.42)
    ) +
    ggplot2::scale_fill_manual(
      name = "Pred1000 > 44",
      values = c("FALSE" = "#2C7BB6", "TRUE" = "#FDE725"),
      labels = c("FALSE", "TRUE")
    ) +
    ggplot2::scale_shape_manual(
      name = "Pred1000 > 44",
      values = c("FALSE" = 21, "TRUE" = 23),
      labels = c("FALSE", "TRUE")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = c(21, 23), color = "black")),
      shape = "none"
    ) +
    ggplot2::labs(x = NULL, y = "Mean rate (0-100 days)") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, color = "black"),
      legend.position = "right",
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )

  save_plot_pair(
    p,
    file.path(figures_dir, output_prefix),
    width = 4.2,
    height = 4.8
  )
  invisible(plot_data)
}

paper_generate_umap_figures <- function(dataset = "invivo",
                                        tables_dir = paper_tables_dir(dataset),
                                        figures_dir = paper_figures_dir(dataset),
                                        support_tables_dir = tables_dir,
                                        figures_wclusters_dir = NULL,
                                        tables_wclusters_dir = NULL,
                                        initial_csv = file.path(tables_dir, paste0(normalize_dataset(dataset), "_deoptim_initial_population.csv")),
                                        best_csv = file.path(tables_dir, paste0(normalize_dataset(dataset), "_best_params_by_seed.csv")),
                                        objective_seed_dir = default_dataset_input_dir(dataset),
                                        output_prefix = paste0(normalize_dataset(dataset), "_deoptim_initial_vs_best_umap"),
                                        best_output_prefix = paste0(normalize_dataset(dataset), "_best_params_umap"),
                                        run_combined = TRUE,
                                        run_sampled = TRUE,
                                        run_best_only = TRUE,
                                        run_clustered_umaps = FALSE,
                                        drop_parameter_table_initial = identical(normalize_dataset(dataset), "invivo"),
                                        shape_by_pred = identical(normalize_dataset(dataset), "invivo"),
                                        umap_seed = 123L,
                                        n_neighbors = 80L,
                                        min_dist = 0.1,
                                        n_threads = max(1L, min(8L, parallel::detectCores(logical = TRUE) %||% 1L)),
                                        pca_n = 10L,
                                        sample_initial_n = 500L,
                                        sample_initial_seed = 123L,
                                        sample_initial_by_seed = TRUE,
                                        cluster_seed = 123L,
                                        cluster_k_min = 2L,
                                        cluster_k_max = 8L,
                                        cluster_silhouette_sample_n = 5000L) {
  dataset <- normalize_dataset(dataset)
  for (pkg in c("ggplot2", "uwot")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Required R package is not installed: ", pkg)
    }
  }

  tables_dir <- normalizePath(path.expand(tables_dir), mustWork = FALSE)
  figures_dir <- normalizePath(path.expand(figures_dir), mustWork = FALSE)
  support_tables_dir <- normalizePath(path.expand(support_tables_dir), mustWork = FALSE)
  initial_csv <- normalizePath(path.expand(initial_csv), mustWork = FALSE)
  best_csv <- normalizePath(path.expand(best_csv), mustWork = FALSE)
  objective_seed_dir <- normalizePath(path.expand(objective_seed_dir), mustWork = FALSE)
  if (isTRUE(run_clustered_umaps)) {
    figures_wclusters_dir <- normalizePath(path.expand(figures_wclusters_dir %||% file.path(dirname(figures_dir), "FiguresWclusters")), mustWork = FALSE)
    tables_wclusters_dir <- normalizePath(path.expand(tables_wclusters_dir %||% file.path(dirname(support_tables_dir), "TablesWclusters")), mustWork = FALSE)
  } else {
    figures_wclusters_dir <- NULL
    tables_wclusters_dir <- NULL
  }
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(support_tables_dir, recursive = TRUE, showWarnings = FALSE)

  if (!run_combined && !run_sampled && !run_best_only) {
    stop("Nothing to run: set run_combined, run_sampled, and/or run_best_only to TRUE.")
  }
  if ((run_combined || run_sampled) && !file.exists(initial_csv)) stop("Missing initial population CSV: ", initial_csv)
  if (!file.exists(best_csv)) stop("Missing best-parameter CSV: ", best_csv)

  params <- umap_parameter_set(dataset)
  log10_params <- umap_log10_parameter_set(dataset)
  expected_n_params <- if (identical(dataset, "invitro")) 14L else 18L
  if (length(params) != expected_n_params) stop("Internal error: UMAP parameter list is not length ", expected_n_params, ".")

  message("Reading best parameters: ", best_csv)
  best_df <- attach_objective(read_csv_plain(best_csv), objective_seed_dir = objective_seed_dir)
  if (!"seed" %in% names(best_df)) stop("Best-parameter CSV must contain a seed column.")
  if (isTRUE(shape_by_pred) && !"pred1000_both_gt44" %in% names(best_df)) {
    stop("Best-parameter CSV must contain pred1000_both_gt44. Regenerate it with util.R first.")
  }
  if (isTRUE(shape_by_pred)) {
    best_df$pred1000_both_gt44 <- coerce_logical_column(best_df$pred1000_both_gt44, "pred1000_both_gt44")
  }
  best_features <- transform_umap_features(best_df, params, log10_params)

  if (run_combined || run_sampled) {
    message("Reading initial population: ", initial_csv)
    initial_df <- read_csv_plain(initial_csv)
    if (!"seed" %in% names(initial_df)) stop("Initial population CSV must contain a seed column.")
    if (isTRUE(drop_parameter_table_initial)) {
      initial_df <- drop_parameter_table_initial_rows(initial_df)
    }
    initial_features <- transform_umap_features(initial_df, params, log10_params)
  }

  output_figure_prefix <- file.path(figures_dir, output_prefix)
  output_table_prefix <- file.path(support_tables_dir, output_prefix)
  pca_figure_prefix <- default_pca_output_prefix(output_figure_prefix)
  pca_table_prefix <- default_pca_output_prefix(output_table_prefix)
  sampled_figure_prefix <- default_sampled_output_prefix(output_figure_prefix, sample_initial_n)
  sampled_table_prefix <- default_sampled_output_prefix(output_table_prefix, sample_initial_n)
  best_figure_prefix <- file.path(figures_dir, best_output_prefix)
  best_table_prefix <- file.path(support_tables_dir, best_output_prefix)
  best_pca_figure_prefix <- default_pca_output_prefix(best_figure_prefix)
  best_pca_table_prefix <- default_pca_output_prefix(best_table_prefix)

  if (run_combined) {
    feature_mat <- standardize_features(rbind(initial_features, best_features))
    direct_emb <- run_umap_embedding(feature_mat, "direct", umap_seed, n_neighbors, min_dist, n_threads)
    direct_plot_data <- build_plot_data(direct_emb, initial_df, best_df, reduction_label = "direct_umap", shape_by_pred = shape_by_pred)
    write_umap_outputs(
      direct_plot_data,
      output_figure_prefix,
      output_table_prefix,
      shape_by_pred = shape_by_pred,
      cluster_feature_mat = feature_mat,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )

    pca_features <- run_pca(
      feature_mat,
      pca_n = pca_n,
      variance_path = paste0(pca_table_prefix, "_pca_variance.csv"),
      variance_figure_prefix = paste0(pca_figure_prefix, "_pca_variance_bar")
    )
    pca_emb <- run_umap_embedding(pca_features, "PCA", umap_seed, n_neighbors, min_dist, n_threads)
    pca_plot_data <- build_plot_data(
      pca_emb,
      initial_df,
      best_df,
      reduction_label = paste0("pca", ncol(pca_features), "_umap"),
      shape_by_pred = shape_by_pred
    )
    write_umap_outputs(
      pca_plot_data,
      pca_figure_prefix,
      pca_table_prefix,
      shape_by_pred = shape_by_pred,
      cluster_feature_mat = pca_features,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )
    message("PCA UMAP retained PCs: ", ncol(pca_features))
  }

  if (run_sampled) {
    sampled_idx <- sample_initial_rows(
      initial_df,
      sample_initial_n,
      sample_initial_seed,
      by_seed = sample_initial_by_seed
    )
    sampled_initial_df <- initial_df[sampled_idx, , drop = FALSE]
    sampled_features <- rbind(initial_features[sampled_idx, , drop = FALSE], best_features)
    sampled_feature_mat <- standardize_features(sampled_features)
    sampled_emb <- run_umap_embedding(
      sampled_feature_mat,
      paste0("sampled-initial-", length(sampled_idx)),
      umap_seed,
      n_neighbors,
      min_dist,
      n_threads
    )
    sampled_plot_data <- build_plot_data(
      sampled_emb,
      sampled_initial_df,
      best_df,
      reduction_label = paste0("sampled", length(sampled_idx), "_umap"),
      shape_by_pred = shape_by_pred
    )
    write_umap_outputs(
      sampled_plot_data,
      sampled_figure_prefix,
      sampled_table_prefix,
      initial_size = 0.6,
      shape_by_pred = shape_by_pred,
      cluster_feature_mat = sampled_feature_mat,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )
    message("Sampled initial UMAP initial points: ", length(sampled_idx))
  }

  if (run_best_only) {
    empty_initial_df <- best_df[0L, , drop = FALSE]
    best_feature_mat <- standardize_features(best_features)
    best_emb <- run_umap_embedding(best_feature_mat, "best-only direct", umap_seed, n_neighbors, min_dist, n_threads)
    best_plot_data <- build_plot_data(best_emb, empty_initial_df, best_df, reduction_label = "best_only_umap", shape_by_pred = shape_by_pred)
    write_umap_outputs(
      best_plot_data,
      best_figure_prefix,
      best_table_prefix,
      best_size = 1.6,
      shape_by_pred = shape_by_pred,
      cluster_feature_mat = best_feature_mat,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )

    best_pca_features <- run_pca(
      best_feature_mat,
      pca_n = pca_n,
      variance_path = paste0(best_pca_table_prefix, "_pca_variance.csv"),
      variance_figure_prefix = paste0(best_pca_figure_prefix, "_pca_variance_bar")
    )
    best_pca_emb <- run_umap_embedding(best_pca_features, "best-only PCA", umap_seed, n_neighbors, min_dist, n_threads)
    best_pca_plot_data <- build_plot_data(
      best_pca_emb,
      empty_initial_df,
      best_df,
      reduction_label = paste0("best_only_pca", ncol(best_pca_features), "_umap"),
      shape_by_pred = shape_by_pred
    )
    write_umap_outputs(
      best_pca_plot_data,
      best_pca_figure_prefix,
      best_pca_table_prefix,
      best_size = 1.6,
      shape_by_pred = shape_by_pred,
      cluster_feature_mat = best_pca_features,
      figures_wclusters_dir = figures_wclusters_dir,
      tables_wclusters_dir = tables_wclusters_dir,
      cluster_seed = cluster_seed,
      cluster_k_min = cluster_k_min,
      cluster_k_max = cluster_k_max,
      cluster_silhouette_sample_n = cluster_silhouette_sample_n
    )
    message("Best-only PCA UMAP retained PCs: ", ncol(best_pca_features))
  }

  message("UMAP parameters: ", paste(params, collapse = ", "))
  invisible(TRUE)
}

read_fit_args_from_seed <- function(seed_dir) {
  path <- file.path(seed_dir, "run_effective_args.tsv")
  if (!file.exists(path)) stop("Missing run_effective_args.tsv: ", path)
  tab <- read_tsv(path)
  req <- c("source", "key", "value")
  if (!all(req %in% names(tab))) stop("Invalid run_effective_args.tsv: ", path)
  tab <- tab[tab$source == "fit_command", , drop = FALSE]
  vals <- as.list(as.character(tab$value))
  names(vals) <- as.character(tab$key)
  vals
}

confirm_natural_to_transformed_initialpop <- function(initial_csv, seed, original_seed_dir) {
  if (!file.exists(initial_csv)) stop("Missing initial population CSV: ", initial_csv)
  pop <- read_csv_plain(initial_csv)
  if (!"seed" %in% names(pop)) stop("Initial population CSV has no seed column: ", initial_csv)
  pop <- pop[as.integer(pop$seed) == as.integer(seed), , drop = FALSE]
  if (!nrow(pop)) stop("No rows found for seed ", seed, " in ", initial_csv)

  param_tab <- read_param_table(original_seed_dir, dataset = "invivo")
  opt_tab <- param_tab[param_tab$estimate, , drop = FALSE]
  transformed <- matrix(NA_real_, nrow = nrow(pop), ncol = nrow(opt_tab))
  colnames(transformed) <- as.character(opt_tab$param_name)
  for (j in seq_len(nrow(opt_tab))) {
    natural_name <- as.character(opt_tab$param_prototype[[j]])
    transformed_name <- as.character(opt_tab$param_name[[j]])
    if (!natural_name %in% names(pop)) {
      stop("Initial population CSV is missing natural parameter column: ", natural_name)
    }
    vals <- suppressWarnings(as.numeric(pop[[natural_name]]))
    if (any(!is.finite(vals))) {
      stop("Initial population has non-finite values for ", natural_name)
    }
    if (startsWith(transformed_name, "log10_")) {
      if (any(vals <= 0)) stop("Cannot log10-transform non-positive values for ", natural_name)
      vals <- log10(vals)
    }
    transformed[, j] <- vals
  }
  transformed
}

confirm_make_comparison_table <- function(original_seed_dir, confirm_seed_dir, output_dir) {
  original_params <- read_best_params(original_seed_dir)
  confirm_params <- read_best_params(confirm_seed_dir)
  param_names <- c(names(original_params), setdiff(names(confirm_params), names(original_params)))

  rows <- vector("list", length(param_names))
  for (i in seq_along(param_names)) {
    nm <- param_names[[i]]
    old <- unname(original_params[[nm]])
    new <- unname(confirm_params[[nm]])
    rows[[i]] <- data.frame(
      parameter_or_metric = nm,
      original = old,
      specified_initial_population = new,
      difference = new - old,
      stringsAsFactors = FALSE
    )
  }

  original_summary <- metric_map(file.path(original_seed_dir, "fit_summary.tsv"))
  confirm_summary <- metric_map(file.path(confirm_seed_dir, "fit_summary.tsv"))
  metric_names <- c(
    "objective",
    "optimizer_deoptim_objective",
    "optimizer_local_objective",
    "optimizer_iter_completed",
    "optimizer_iter_target",
    "optimizer_local_accepted",
    "optimizer_local_convergence"
  )
  metric_rows <- lapply(metric_names, function(nm) {
    old_raw <- metric_value(original_summary, nm)
    new_raw <- metric_value(confirm_summary, nm)
    old_num <- suppressWarnings(as.numeric(old_raw))
    new_num <- suppressWarnings(as.numeric(new_raw))
    numeric_pair <- is.finite(old_num) && is.finite(new_num)
    data.frame(
      parameter_or_metric = nm,
      original = if (numeric_pair) old_num else old_raw,
      specified_initial_population = if (numeric_pair) new_num else new_raw,
      difference = if (numeric_pair) new_num - old_num else if (identical(old_raw, new_raw)) "0" else NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, c(rows, metric_rows))
  out_path <- file.path(output_dir, "seed1_initial_population_reproducibility_comparison.csv")
  utils::write.csv(out, out_path, quote = FALSE, row.names = FALSE)
  message("Wrote comparison table: ", out_path)
  invisible(out)
}

confirm_invivo_seed_initial_population_fit <- function(seed = 1L,
                                                       original_seed_dir = file.path(DEFAULT_INVIVO_INPUT_DIR, paste0("seed", as.integer(seed))),
                                                       initial_population_csv = file.path(paper_tables_dir("invivo"), "invivo_deoptim_initial_population.csv"),
                                                       output_dir = file.path(paper_tables_dir("invivo"), "ConfirmFits"),
                                                       backend = NULL) {
  seed <- as.integer(seed)
  original_seed_dir <- normalizePath(path.expand(original_seed_dir), mustWork = FALSE)
  initial_population_csv <- normalizePath(path.expand(initial_population_csv), mustWork = FALSE)
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  confirm_seed_dir <- file.path(output_dir, paste0("seed", seed, "_from_initial_population"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(confirm_seed_dir)) unlink(confirm_seed_dir, recursive = TRUE, force = TRUE)
  dir.create(confirm_seed_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(original_seed_dir)) stop("Original seed directory does not exist: ", original_seed_dir)

  this_script_dir <- script_dir()
  workflow_root <- normalizePath(file.path(dirname(this_script_dir), "O2_supply_demand_MAP"), mustWork = FALSE)
  optimizer_dir <- file.path(workflow_root, "optimizer")
  backend_path <- normalizePath(
    backend %||% file.path(workflow_root, "util", "o2_supply_demand_map_fit_invivo_backend.R"),
    mustWork = FALSE
  )
  if (!file.exists(backend_path)) stop("Backend script not found: ", backend_path)

  initialpop_t <- confirm_natural_to_transformed_initialpop(
    initial_csv = initial_population_csv,
    seed = seed,
    original_seed_dir = original_seed_dir
  )

  fit_args <- read_fit_args_from_seed(original_seed_dir)
  fit_args$seed <- as.character(seed)
  fit_args$out_dir <- confirm_seed_dir
  fit_args$mode <- "fit_seed"
  fit_args$fit_invivo <- "TRUE"

  raw_command_args <- base::commandArgs
  backend_env <- new.env(parent = globalenv())
  backend_env$commandArgs <- function(trailingOnly = FALSE) {
    args <- raw_command_args(trailingOnly = trailingOnly)
    if (!isTRUE(trailingOnly)) {
      args <- args[!grepl("^--file=", args)]
    }
    args
  }
  sys.source(backend_path, envir = backend_env, chdir = TRUE)
  backend_get_script_dir <- function(default = getwd()) {
    normalizePath(optimizer_dir, mustWork = FALSE)
  }
  assign("get_script_dir", backend_get_script_dir, envir = backend_env)
  assign("o2sd_get_script_dir", backend_get_script_dir, envir = backend_env)
  environment(backend_env$main_fit_single_seed)$get_script_dir <- backend_get_script_dir
  environment(backend_env$main)$get_script_dir <- backend_get_script_dir
  message("Backend path: ", backend_path)
  message("Backend optimizer dir: ", backend_get_script_dir())

  original_build_de_initialpop <- get(".build_de_initialpop", envir = backend_env, inherits = FALSE)
  assign(
    ".build_de_initialpop",
    function(np, lower, upper, init_use = NULL, mode = "hybrid", uniform_frac = 0.3, sigma_frac = 0.1) {
      regenerated_plan <- original_build_de_initialpop(
        np = np,
        lower = lower,
        upper = upper,
        init_use = init_use,
        mode = mode,
        uniform_frac = uniform_frac,
        sigma_frac = sigma_frac
      )
      regenerated <- if (is.list(regenerated_plan) && !is.null(regenerated_plan$pop)) {
        regenerated_plan$pop
      } else {
        regenerated_plan
      }
      requested <- initialpop_t[, names(lower), drop = FALSE]
      if (nrow(requested) != as.integer(np)) {
        stop("Initial population row count mismatch: CSV has ", nrow(requested), ", NP is ", np)
      }
      if (any(!is.finite(requested))) stop("Initial population contains non-finite transformed values")
      if (any(requested < matrix(lower, nrow = nrow(requested), ncol = length(lower), byrow = TRUE) - 1e-8) ||
          any(requested > matrix(upper, nrow = nrow(requested), ncol = length(upper), byrow = TRUE) + 1e-8)) {
        stop("Initial population contains values outside transformed optimizer bounds.")
      }
      if (is.matrix(regenerated) && all(names(lower) %in% colnames(regenerated))) {
        max_abs_delta <- max(abs(requested - regenerated[, names(lower), drop = FALSE]))
        message("Using supplied initial population; max abs delta vs regenerated population = ", signif(max_abs_delta, 17))
      } else {
        message("Using supplied initial population; regenerated population was not matrix-like for delta check.")
      }
      if (is.list(regenerated_plan)) {
        regenerated_plan$pop <- requested
        regenerated_plan$mode_effective <- paste0("specified_", regenerated_plan$mode_effective %||% "initial_population")
        return(regenerated_plan)
      }
      list(
        pop = requested,
        mode_effective = "specified_initial_population",
        warm_start_used = !is.null(init_use),
        n_local = NA_integer_,
        n_uniform = NA_integer_
      )
    },
    envir = backend_env
  )

  message("Running confirm fit into: ", confirm_seed_dir)
  backend_env$main_fit_single_seed(fit_args)
  comp <- confirm_make_comparison_table(
    original_seed_dir = original_seed_dir,
    confirm_seed_dir = confirm_seed_dir,
    output_dir = output_dir
  )
  diff_num <- suppressWarnings(as.numeric(comp$difference))
  finite_diff <- diff_num[is.finite(diff_num)]
  if (length(finite_diff)) {
    message("Max absolute numeric difference: ", signif(max(abs(finite_diff)), 17))
  }
  invisible(comp)
}
