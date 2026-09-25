#!/usr/bin/env Rscript

.script_dir <- local({
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
  if (length(frame_files) > 0L) dirname(frame_files[[length(frame_files)]]) else getwd()
})
SCRIPT_DIR <- normalizePath(.script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
rm(.script_dir)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_fit_results_utils.R"), local = TRUE)
`%||%` <- o2fr_null_coalesce

TARGET_PARAMETERS <- c("p_misseg", "O2_crit", "alpha_o2", "mu_hp")
CONTEXT_LEVELS <- c("in vivo", "in vitro")
OBJECTIVE_TYPES <- c("objective", "objective_invivo", "objective_invitro")
OUTPUT_PREFIX <- "joint_sigma_soft_coupled_paired_seed_comparison"

parse_args <- o2fr_parse_args
as_num <- o2fr_as_num

as_bool <- o2fr_as_bool

 seed_order_key <- function(x) {
  x <- basename(as.character(x))
  m <- regexec("^seed([0-9]+)$", x)
  hit <- regmatches(x, m)
  out <- rep(Inf, length(x))
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) == 2L) out[[i]] <- suppressWarnings(as.numeric(hit[[i]][[2]]))
  }
  out
}

format_sigma <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || !is.finite(x[[1]])) return("unknown")
  x <- x[[1]]
  if (abs(x) >= 1000 || (abs(x) > 0 && abs(x) < 0.001)) {
    format(x, scientific = TRUE, digits = 4, trim = TRUE)
  } else {
    format(x, scientific = FALSE, digits = 4, trim = TRUE)
  }
}

read_metric_map_optional <- function(path) {
  if (!file.exists(path)) return(character(0))
  tab <- tryCatch(
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !all(c("metric", "value") %in% names(tab))) return(character(0))
  setNames(tab$value, tab$metric)
}

resolve_default_results_dir <- function() {
  normalizePath(file.path(WORKFLOW_ROOT, "..", "..", "results"), mustWork = FALSE)
}

resolve_run_dirs <- function(results_dir, run_dirs_arg = NULL) {
  if (!is.null(run_dirs_arg) && nzchar(trimws(run_dirs_arg))) {
    raw <- trimws(unlist(strsplit(run_dirs_arg, ",", fixed = TRUE)))
    raw <- raw[nzchar(raw)]
    dirs <- vapply(raw, function(x) {
      candidate <- x
      if (!dir.exists(candidate)) candidate <- file.path(results_dir, x)
      normalizePath(candidate, mustWork = TRUE)
    }, character(1))
    return(unname(dirs))
  }

  dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
  keep <- vapply(dirs, function(d) {
    seed_dirs <- list.dirs(d, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    length(seed_dirs) > 0L && any(file.exists(file.path(seed_dirs, "joint_soft_coupling.tsv")))
  }, logical(1))
  dirs <- dirs[keep]
  if (!length(dirs)) {
    stop("No run directories with seed*/joint_soft_coupling.tsv found under: ", results_dir, call. = FALSE)
  }
  normalizePath(dirs[order(basename(dirs))], mustWork = TRUE)
}

read_seed_soft_table <- function(run_dir, seed_dir, target_parameters = TARGET_PARAMETERS) {
  soft_path <- file.path(seed_dir, "joint_soft_coupling.tsv")
  if (!file.exists(soft_path)) return(NULL)
  soft_tab <- tryCatch(
    utils::read.delim(soft_path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(soft_tab) || !nrow(soft_tab) || !"parameter" %in% names(soft_tab)) return(NULL)
  soft_tab <- soft_tab[as.character(soft_tab$parameter) %in% target_parameters, , drop = FALSE]
  if (!nrow(soft_tab)) return(NULL)

  seed <- basename(seed_dir)
  fit_summary <- read_metric_map_optional(file.path(seed_dir, "fit_summary.tsv"))
  run_sigma <- as_num(fit_summary[["joint_soft_coupling_sigma_default"]], NA_real_)
  if (!is.finite(run_sigma) && "regularization_sigma" %in% names(soft_tab)) {
    run_sigma <- as_num(soft_tab$regularization_sigma[is.finite(suppressWarnings(as.numeric(soft_tab$regularization_sigma)))][[1]], NA_real_)
  }

  get_num <- function(col) {
    if (col %in% names(soft_tab)) suppressWarnings(as.numeric(soft_tab[[col]])) else rep(NA_real_, nrow(soft_tab))
  }
  get_chr <- function(col) {
    if (col %in% names(soft_tab)) as.character(soft_tab[[col]]) else rep(NA_character_, nrow(soft_tab))
  }
  get_bool <- function(col) {
    if (!(col %in% names(soft_tab))) return(rep(FALSE, nrow(soft_tab)))
    vapply(soft_tab[[col]], as_bool, logical(1), default = FALSE)
  }

  vivo_natural <- get_num("vivo_natural")
  vitro_natural <- get_num("vitro_natural")
  ratio <- get_num("ratio_vivo_to_vitro")
  missing_ratio <- !is.finite(ratio) & is.finite(vivo_natural) & is.finite(vitro_natural) & vitro_natural != 0
  ratio[missing_ratio] <- vivo_natural[missing_ratio] / vitro_natural[missing_ratio]

  data.frame(
    run_label = basename(run_dir),
    run_dir = normalizePath(run_dir, mustWork = FALSE),
    seed = seed,
    seed_id = seed_order_key(seed),
    parameter = as.character(soft_tab$parameter),
    parameter_order = match(as.character(soft_tab$parameter), target_parameters),
    sigma = run_sigma,
    sigma_label = paste0("sigma=", format_sigma(run_sigma)),
    regularization_sigma = get_num("regularization_sigma"),
    penalty_type = get_chr("penalty_type"),
    welsch_c = get_num("welsch_c"),
    welsch_transition_delta = get_num("welsch_transition_delta"),
    abs_delta_over_sigma = get_num("abs_delta_over_sigma"),
    welsch_saturation_fraction = get_num("welsch_saturation_fraction"),
    penalty_region = get_chr("penalty_region"),
    split_enabled = get_bool("split_enabled"),
    center_name = get_chr("center_name"),
    delta_name = get_chr("delta_name"),
    center_transformed = get_num("center_transformed"),
    delta_transformed = get_num("delta_transformed"),
    vivo_transformed = get_num("vivo_transformed"),
    vitro_transformed = get_num("vitro_transformed"),
    center_natural = get_num("center_natural"),
    vivo_natural = vivo_natural,
    vitro_natural = vitro_natural,
    delta_interpretable = get_num("delta_interpretable"),
    natural_difference_vivo_to_vitro = get_num("natural_difference_vivo_to_vitro"),
    transformed_difference_vivo_to_vitro = get_num("transformed_difference_vivo_to_vitro"),
    log10_ratio_vivo_to_vitro = get_num("log10_ratio_vivo_to_vitro"),
    fold_change_vivo_to_vitro = get_num("fold_change_vivo_to_vitro"),
    ratio_vivo_to_vitro = ratio,
    odds_ratio_vivo_to_vitro = get_num("odds_ratio_vivo_to_vitro"),
    penalty_paid = get_num("penalty_paid"),
    joint_union_lower_transformed = get_num("joint_union_lower_transformed"),
    joint_union_upper_transformed = get_num("joint_union_upper_transformed"),
    feasible_at_solution = get_bool("feasible_at_solution"),
    objective = as_num(fit_summary[["objective"]], NA_real_),
    objective_invivo = as_num(fit_summary[["objective_invivo"]], NA_real_),
    objective_invitro = as_num(fit_summary[["objective_invitro"]], NA_real_),
    objective_soft_coupling = as_num(fit_summary[["objective_soft_coupling"]], NA_real_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_all_soft_tables <- function(run_dirs, target_parameters = TARGET_PARAMETERS) {
  rows <- list()
  k <- 0L
  for (run_dir in run_dirs) {
    seed_dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    seed_dirs <- seed_dirs[order(seed_order_key(basename(seed_dirs)), basename(seed_dirs))]
    for (seed_dir in seed_dirs) {
      row <- read_seed_soft_table(run_dir, seed_dir, target_parameters = target_parameters)
      if (is.null(row) || !nrow(row)) next
      k <- k + 1L
      rows[[k]] <- row
    }
  }
  if (!length(rows)) {
    stop("No soft-coupling rows were read from run directories.", call. = FALSE)
  }
  out <- do.call(rbind, rows)
  out$parameter <- factor(as.character(out$parameter), levels = target_parameters)
  out <- out[order(out$parameter_order, out$seed_id, out$seed, out$sigma, out$run_label), , drop = FALSE]
  row.names(out) <- NULL
  out
}

filter_complete_paired_seeds <- function(df, target_parameters = TARGET_PARAMETERS) {
  run_labels <- unique(as.character(df$run_label))
  needed <- length(run_labels) * length(target_parameters)
  keys <- unique(df[, c("seed", "run_label", "parameter"), drop = FALSE])
  counts <- stats::aggregate(
    run_label ~ seed,
    data = keys,
    FUN = length
  )
  names(counts)[names(counts) == "run_label"] <- "complete_cells"
  complete_seeds <- as.character(counts$seed[counts$complete_cells == needed])
  df <- df[as.character(df$seed) %in% complete_seeds, , drop = FALSE]
  df <- df[order(df$parameter_order, df$seed_id, df$seed, df$sigma, df$run_label), , drop = FALSE]
  row.names(df) <- NULL
  attr(df, "complete_seed_count") <- length(unique(as.character(df$seed)))
  attr(df, "all_seed_count") <- length(unique(as.character(keys$seed)))
  df
}

sigma_levels_from_df <- function(df) {
  sigma_map <- unique(df[, c("sigma", "sigma_label"), drop = FALSE])
  sigma_map <- sigma_map[order(sigma_map$sigma, sigma_map$sigma_label, na.last = TRUE), , drop = FALSE]
  as.character(sigma_map$sigma_label)
}

make_value_long <- function(df, sigma_levels) {
  vivo <- data.frame(
    df[, c("run_label", "run_dir", "seed", "seed_id", "parameter", "parameter_order", "sigma", "sigma_label"), drop = FALSE],
    context = "in vivo",
    value = df$vivo_natural,
    stringsAsFactors = FALSE
  )
  vitro <- data.frame(
    df[, c("run_label", "run_dir", "seed", "seed_id", "parameter", "parameter_order", "sigma", "sigma_label"), drop = FALSE],
    context = "in vitro",
    value = df$vitro_natural,
    stringsAsFactors = FALSE
  )
  out <- rbind(vivo, vitro)
  out$context <- factor(out$context, levels = CONTEXT_LEVELS)
  out$sigma_label <- factor(as.character(out$sigma_label), levels = sigma_levels)
  group_levels <- as.vector(t(outer(CONTEXT_LEVELS, sigma_levels, paste, sep = "__")))
  out$value_group <- factor(
    paste(as.character(out$context), as.character(out$sigma_label), sep = "__"),
    levels = group_levels
  )
  out
}

make_objective_long <- function(df, sigma_levels) {
  seed_df <- unique(df[, c(
    "run_label",
    "run_dir",
    "seed",
    "seed_id",
    "sigma",
    "sigma_label",
    OBJECTIVE_TYPES
  ), drop = FALSE])
  seed_df <- seed_df[order(seed_df$seed_id, seed_df$seed, seed_df$sigma), , drop = FALSE]
  rows <- lapply(OBJECTIVE_TYPES, function(objective_type) {
    data.frame(
      seed_df[, c("run_label", "run_dir", "seed", "seed_id", "sigma", "sigma_label"), drop = FALSE],
      objective_type = objective_type,
      objective_value = suppressWarnings(as.numeric(seed_df[[objective_type]])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$objective_type <- factor(out$objective_type, levels = OBJECTIVE_TYPES)
  out$sigma_label <- factor(as.character(out$sigma_label), levels = sigma_levels)
  group_levels <- as.vector(t(outer(OBJECTIVE_TYPES, sigma_levels, paste, sep = "__")))
  out$objective_group <- factor(
    paste(as.character(out$objective_type), as.character(out$sigma_label), sep = "__"),
    levels = group_levels
  )
  out <- out[order(out$objective_type, out$seed_id, out$seed, out$sigma), , drop = FALSE]
  row.names(out) <- NULL
  out
}

 build_two_sigma_contrast <- function(df) {
  sigma_levels <- sigma_levels_from_df(df)
  if (length(sigma_levels) != 2L) return(data.frame())
  low_label <- sigma_levels[[1]]
  high_label <- sigma_levels[[2]]
  low <- df[as.character(df$sigma_label) == low_label, , drop = FALSE]
  high <- df[as.character(df$sigma_label) == high_label, , drop = FALSE]
  key <- c("seed", "parameter")
  merged <- merge(
    low,
    high,
    by = key,
    suffixes = c("_low_sigma", "_high_sigma"),
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(merged)) return(data.frame())
  ratio_safe <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)
  out <- data.frame(
    seed = merged$seed,
    seed_id = merged$seed_id_low_sigma,
    parameter = merged$parameter,
    parameter_order = merged$parameter_order_low_sigma,
    sigma_low = merged$sigma_low_sigma,
    sigma_high = merged$sigma_high_sigma,
    sigma_label_low = merged$sigma_label_low_sigma,
    sigma_label_high = merged$sigma_label_high_sigma,
    run_low = merged$run_label_low_sigma,
    run_high = merged$run_label_high_sigma,
    vivo_natural_low_sigma = merged$vivo_natural_low_sigma,
    vivo_natural_high_sigma = merged$vivo_natural_high_sigma,
    vivo_natural_delta_high_minus_low = merged$vivo_natural_high_sigma - merged$vivo_natural_low_sigma,
    vivo_natural_ratio_high_to_low = ratio_safe(merged$vivo_natural_high_sigma, merged$vivo_natural_low_sigma),
    vitro_natural_low_sigma = merged$vitro_natural_low_sigma,
    vitro_natural_high_sigma = merged$vitro_natural_high_sigma,
    vitro_natural_delta_high_minus_low = merged$vitro_natural_high_sigma - merged$vitro_natural_low_sigma,
    vitro_natural_ratio_high_to_low = ratio_safe(merged$vitro_natural_high_sigma, merged$vitro_natural_low_sigma),
    ratio_vivo_to_vitro_low_sigma = merged$ratio_vivo_to_vitro_low_sigma,
    ratio_vivo_to_vitro_high_sigma = merged$ratio_vivo_to_vitro_high_sigma,
    ratio_vivo_to_vitro_delta_high_minus_low = merged$ratio_vivo_to_vitro_high_sigma - merged$ratio_vivo_to_vitro_low_sigma,
    ratio_vivo_to_vitro_ratio_high_to_low = ratio_safe(merged$ratio_vivo_to_vitro_high_sigma, merged$ratio_vivo_to_vitro_low_sigma),
    feasible_low_sigma = merged$feasible_at_solution_low_sigma,
    feasible_high_sigma = merged$feasible_at_solution_high_sigma,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[order(out$parameter_order, out$seed_id, out$seed), , drop = FALSE]
  row.names(out) <- NULL
  out
}

 main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  results_dir <- normalizePath(args$results_dir %||% resolve_default_results_dir(), mustWork = TRUE)
  out_dir <- normalizePath(args$out_dir %||% results_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_dirs <- resolve_run_dirs(results_dir, args$run_dirs %||% NULL)
  if (length(run_dirs) < 2L) {
    stop("Need at least two run directories for paired sigma comparison.", call. = FALSE)
  }

  message("Reading run directories: ", paste(basename(run_dirs), collapse = ", "))
  raw_df <- read_all_soft_tables(run_dirs, target_parameters = TARGET_PARAMETERS)
  paired_df <- filter_complete_paired_seeds(raw_df, target_parameters = TARGET_PARAMETERS)
  if (!nrow(paired_df)) {
    stop("No complete paired seeds remained after requiring all target parameters in all runs.", call. = FALSE)
  }
  all_seed_count <- attr(paired_df, "all_seed_count")
  complete_seed_count <- attr(paired_df, "complete_seed_count")
  if (!is.null(all_seed_count) && !is.null(complete_seed_count) && all_seed_count != complete_seed_count) {
    message("Keeping ", complete_seed_count, " complete paired seeds out of ", all_seed_count, " seeds observed in any run.")
  } else {
    message("Keeping ", length(unique(as.character(paired_df$seed))), " complete paired seeds.")
  }

  sigma_levels <- sigma_levels_from_df(paired_df)
  paired_df$sigma_label <- factor(as.character(paired_df$sigma_label), levels = sigma_levels)
  paired_df$parameter <- factor(as.character(paired_df$parameter), levels = TARGET_PARAMETERS)
  paired_df <- paired_df[order(paired_df$parameter_order, paired_df$seed_id, paired_df$seed, paired_df$sigma), , drop = FALSE]

  summary_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_summary.tsv"))
  utils::write.table(paired_df, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  contrast_df <- build_two_sigma_contrast(paired_df)
  contrast_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_two_sigma_contrast.tsv"))
  if (nrow(contrast_df)) {
    utils::write.table(contrast_df, file = contrast_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  value_long <- make_value_long(paired_df, sigma_levels = sigma_levels)
  value_long_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_value_long.tsv"))
  utils::write.table(value_long, file = value_long_path, sep = "\t", quote = FALSE, row.names = FALSE)

  objective_long <- make_objective_long(paired_df, sigma_levels = sigma_levels)
  objective_long_path <- file.path(out_dir, paste0(OUTPUT_PREFIX, "_objective_long.tsv"))
  utils::write.table(objective_long, file = objective_long_path, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(stage = "analysis", file = basename(c(summary_path, contrast_path, value_long_path, objective_long_path)), stringsAsFactors = FALSE), file.path(out_dir, "analysis_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote summary table: ", summary_path)
  if (nrow(contrast_df)) message("Wrote paired contrast table: ", contrast_path)
  message("Wrote value-long table: ", value_long_path)
  message("Wrote objective-long table: ", objective_long_path)
}

if (sys.nframe() == 0L) {
  main()
}
