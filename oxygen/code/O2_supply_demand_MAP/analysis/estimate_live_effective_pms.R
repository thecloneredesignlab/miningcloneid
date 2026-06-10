#!/usr/bin/env Rscript

# Usage:
#   Rscript oxygen/code/O2_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
#     --seed_dir=/abs/path/to/seed7
#
#   Rscript oxygen/code/O2_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
#     --run_dir=/abs/path/to/run_dir \
#     --seed=7
#
#   Rscript oxygen/code/O2_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
#     --run_dir=/abs/path/to/run_dir
#
# Optional:
#   --out_dir=/abs/path/to/output_dir
#   --seed_id=seed7

.o2sd_live_pms_bootstrap_script_dir <- local({
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
SCRIPT_DIR <- normalizePath(.o2sd_live_pms_bootstrap_script_dir, mustWork = FALSE)
WORKFLOW_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = FALSE)
source(file.path(WORKFLOW_ROOT, "util", "o2_supply_demand_map_shared.R"), local = environment())
rm(.o2sd_live_pms_bootstrap_script_dir)

`%||%` <- o2sd_null_coalesce
parse_args <- o2sd_parse_args
as_num <- o2sd_as_num
as_int <- o2sd_as_int
as_bool <- o2sd_as_bool

resolve_path_value <- function(path_value, base_dir = getwd()) {
  txt <- path_value
  if (is.null(txt) || !length(txt)) return(NULL)
  txt <- as.character(txt[[1]])
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  if (startsWith(txt, "~")) return(normalizePath(path.expand(txt), mustWork = FALSE))
  if (grepl("^(/|[A-Za-z]:[/\\\\])", txt)) return(normalizePath(txt, mustWork = FALSE))
  normalizePath(file.path(base_dir, txt), mustWork = FALSE)
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file was not found: ", path)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

read_best_params_map <- function(seed_dir) {
  path <- file.path(seed_dir, "best_params.tsv")
  tab <- read_required_tsv(path)
  if (!all(c("parameter", "value") %in% names(tab))) {
    stop("best_params.tsv is missing required columns in ", path)
  }
  vals <- suppressWarnings(as.numeric(tab$value))
  if (any(!is.finite(vals))) {
    stop("Non-finite best_params values were found in ", path)
  }
  setNames(vals, tab$parameter)
}

read_fit_config <- function(seed_dir) {
  path <- file.path(seed_dir, "fit_config.rds")
  if (!file.exists(path)) {
    stop("fit_config.rds was not found: ", path)
  }
  readRDS(path)
}

choose_rank_column <- function(tab) {
  candidates <- c(
    "recommend_rank_burden_ploidy_boundary_first",
    "recommend_rank_ploidy_burden_boundary_first",
    "recommend_rank_ploidy_boundary_first",
    "recommend_rank_ploidy_first",
    "boundary_rank_active_support",
    "objective_rank"
  )
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit) > 0L) return(hit[[1]])
  NULL
}

seed_numeric_key <- function(seed_ids) {
  out <- rep(Inf, length(seed_ids))
  m <- regexec("^seed([0-9]+)$", seed_ids)
  hit <- regmatches(seed_ids, m)
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) == 2L) out[[i]] <- suppressWarnings(as.numeric(hit[[i]][[2]]))
  }
  out
}

resolve_seed_dir <- function(argv) {
  seed_dir <- resolve_path_value(argv$seed_dir, getwd())
  if (!is.null(seed_dir)) {
    seed_dir <- normalizePath(seed_dir, mustWork = TRUE)
    return(list(seed_dir = seed_dir, seed_id = basename(seed_dir), selected_by = "seed_dir"))
  }

  run_dir <- resolve_path_value(argv$run_dir, getwd())
  if (is.null(run_dir)) {
    stop("Either --seed_dir or --run_dir must be supplied.")
  }
  run_dir <- normalizePath(run_dir, mustWork = TRUE)

  seed_id <- NULL
  if (!is.null(argv$seed_id) && nzchar(trimws(as.character(argv$seed_id)))) {
    seed_id <- trimws(as.character(argv$seed_id))
  } else if (!is.null(argv$seed) && nzchar(trimws(as.character(argv$seed)))) {
    seed_num <- trimws(as.character(argv$seed))
    seed_id <- if (startsWith(seed_num, "seed")) seed_num else paste0("seed", seed_num)
  }

  if (!is.null(seed_id)) {
    seed_dir <- file.path(run_dir, seed_id)
    if (!dir.exists(seed_dir)) {
      stop("Requested seed directory was not found: ", seed_dir)
    }
    return(list(seed_dir = normalizePath(seed_dir, mustWork = TRUE), seed_id = seed_id, selected_by = "explicit_seed"))
  }

  summary_path <- file.path(run_dir, "extra_results", "seed_summary.tsv")
  if (!file.exists(summary_path)) {
    stop("No --seed/--seed_id was supplied and extra_results/seed_summary.tsv was not found in ", run_dir)
  }
  tab <- read_required_tsv(summary_path)
  if (!("seed" %in% names(tab))) {
    stop("seed_summary.tsv is missing required column 'seed': ", summary_path)
  }
  rank_col <- choose_rank_column(tab)
  if (!is.null(rank_col)) {
    tab$.__rank__ <- suppressWarnings(as.numeric(tab[[rank_col]]))
    ord <- order(tab$.__rank__, suppressWarnings(as.numeric(tab$objective)), seed_numeric_key(tab$seed), na.last = TRUE)
    tab <- tab[ord, , drop = FALSE]
    keep <- which(is.finite(tab$.__rank__))
    if (!length(keep)) {
      stop("No finite values were found in rank column ", rank_col, " in ", summary_path)
    }
    seed_id <- as.character(tab$seed[[keep[[1]]]])
    selected_by <- paste0("extra_results_rank:", rank_col)
  } else {
    if (!("objective" %in% names(tab))) {
      stop("No rank column or objective column was found in ", summary_path)
    }
    ord <- order(suppressWarnings(as.numeric(tab$objective)), seed_numeric_key(tab$seed), na.last = TRUE)
    tab <- tab[ord, , drop = FALSE]
    keep <- which(is.finite(suppressWarnings(as.numeric(tab$objective))))
    if (!length(keep)) {
      stop("No finite objective values were found in ", summary_path)
    }
    seed_id <- as.character(tab$seed[[keep[[1]]]])
    selected_by <- "extra_results_objective"
  }

  seed_dir <- file.path(run_dir, seed_id)
  if (!dir.exists(seed_dir)) {
    stop("Selected seed directory was not found: ", seed_dir)
  }
  list(seed_dir = normalizePath(seed_dir, mustWork = TRUE), seed_id = seed_id, selected_by = selected_by)
}

compute_h_o2 <- function(O2, O2_crit, n_O) {
  O2_crit_use <- if (is.finite(O2_crit) && O2_crit > 0) O2_crit else 1e-12
  n_O_use <- if (is.finite(n_O) && n_O > 0) n_O else 1.0
  num <- O2_crit_use^n_O_use
  den <- num + (pmax(O2, 0)^n_O_use)
  out <- num / pmax(den, 1e-12)
  pmax(0, pmin(1, out))
}

compute_mu_eff <- function(O2, N, mu_hp, gamma_mu, O2_crit, n_O, N_dip, ploidy_O2_death) {
  mu_hp_use <- if (is.finite(mu_hp) && mu_hp > 0) mu_hp else 0.0
  gamma_mu_use <- if (is.finite(gamma_mu) && gamma_mu > 0) gamma_mu else 1.0
  h_o2 <- compute_h_o2(O2 = O2, O2_crit = O2_crit, n_O = n_O)
  if (identical(ploidy_O2_death, "uniform")) {
    return(mu_hp_use * h_o2)
  }
  if (identical(ploidy_O2_death, "diploid_NULL")) {
    return(mu_hp_use * h_o2 * (1 + pmax(N / N_dip - 1, 0)^gamma_mu_use))
  }
  mu_hp_use * h_o2 * pmax(N / N_dip, 0)^gamma_mu_use
}

compute_p_ms <- function(mu_eff, p_mis_base, p_misseg, k_o_mis) {
  p_base <- if (is.finite(p_mis_base) && p_mis_base >= 0) p_mis_base else 1e-5
  p_amp <- if (is.finite(p_misseg) && p_misseg >= 0) p_misseg else 0.0
  k_use <- if (is.finite(k_o_mis) && k_o_mis > 0) k_o_mis else 1e-12
  frac <- pmax(mu_eff, 0) / pmax(mu_eff + k_use, 1e-12)
  pmax(0, pmin(1, p_base + p_amp * frac))
}

aggregate_sample_day <- function(merged, p_misseg_parameter, N_UNIT) {
  out <- stats::aggregate(
    cbind(
      weighted_fraction = merged$fraction,
      weighted_mean_N = merged$fraction * merged$N,
      weighted_o2 = merged$fraction * merged$o2_pct,
      weighted_mu_eff = merged$fraction * merged$mu_eff,
      weighted_p_ms = merged$fraction * merged$p_ms,
      weighted_p_ms_retained = merged$fraction * merged$p_ms_retained_proxy
    ),
    by = merged[c("harvest", "cohort", "dose", "day")],
    FUN = sum
  )
  out$sample_id <- paste(out$harvest, out$cohort, format(out$dose, trim = TRUE, scientific = FALSE), sep = "__")
  out$live_weighted_mean_N <- out$weighted_mean_N / pmax(out$weighted_fraction, 1e-12)
  out$live_weighted_mean_ploidy <- out$live_weighted_mean_N / pmax(as.numeric(N_UNIT), 1e-12)
  out$live_weighted_o2_pct <- out$weighted_o2 / pmax(out$weighted_fraction, 1e-12)
  out$live_weighted_mu_eff <- out$weighted_mu_eff / pmax(out$weighted_fraction, 1e-12)
  out$live_weighted_effective_p_ms <- out$weighted_p_ms / pmax(out$weighted_fraction, 1e-12)
  out$live_weighted_retained_p_ms_proxy <- out$weighted_p_ms_retained / pmax(out$weighted_fraction, 1e-12)
  out$p_misseg_parameter <- p_misseg_parameter
  out$abs_diff_vs_p_misseg <- p_misseg_parameter - out$live_weighted_effective_p_ms
  out$ratio_vs_p_misseg <- out$live_weighted_effective_p_ms / pmax(p_misseg_parameter, 1e-12)
  out$retained_proxy_abs_diff_vs_p_misseg <- p_misseg_parameter - out$live_weighted_retained_p_ms_proxy
  out$retained_proxy_ratio_vs_p_misseg <- out$live_weighted_retained_p_ms_proxy / pmax(p_misseg_parameter, 1e-12)
  out[order(out$cohort, out$harvest, out$dose, out$day), , drop = FALSE]
}

summarise_sample_day_table <- function(tab, p_misseg_parameter, label) {
  if (!nrow(tab)) {
    return(data.frame(
      summary_scope = label,
      n_points = 0L,
      p_misseg_parameter = p_misseg_parameter,
      live_weighted_effective_p_ms_mean = NA_real_,
      live_weighted_effective_p_ms_median = NA_real_,
      live_weighted_effective_p_ms_min = NA_real_,
      live_weighted_effective_p_ms_max = NA_real_,
      live_weighted_effective_p_ms_sd = NA_real_,
      abs_diff_vs_p_misseg_mean = NA_real_,
      ratio_vs_p_misseg_mean = NA_real_,
      live_weighted_retained_p_ms_proxy_mean = NA_real_,
      live_weighted_retained_p_ms_proxy_median = NA_real_,
      live_weighted_retained_p_ms_proxy_min = NA_real_,
      live_weighted_retained_p_ms_proxy_max = NA_real_,
      retained_proxy_abs_diff_vs_p_misseg_mean = NA_real_,
      retained_proxy_ratio_vs_p_misseg_mean = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  x <- as.numeric(tab$live_weighted_effective_p_ms)
  y <- as.numeric(tab$live_weighted_retained_p_ms_proxy)
  data.frame(
    summary_scope = label,
    n_points = nrow(tab),
    p_misseg_parameter = p_misseg_parameter,
    live_weighted_effective_p_ms_mean = mean(x),
    live_weighted_effective_p_ms_median = stats::median(x),
    live_weighted_effective_p_ms_min = min(x),
    live_weighted_effective_p_ms_max = max(x),
    live_weighted_effective_p_ms_sd = stats::sd(x),
    abs_diff_vs_p_misseg_mean = mean(p_misseg_parameter - x),
    ratio_vs_p_misseg_mean = mean(x) / pmax(p_misseg_parameter, 1e-12),
    live_weighted_retained_p_ms_proxy_mean = mean(y),
    live_weighted_retained_p_ms_proxy_median = stats::median(y),
    live_weighted_retained_p_ms_proxy_min = min(y),
    live_weighted_retained_p_ms_proxy_max = max(y),
    retained_proxy_abs_diff_vs_p_misseg_mean = mean(p_misseg_parameter - y),
    retained_proxy_ratio_vs_p_misseg_mean = mean(y) / pmax(p_misseg_parameter, 1e-12),
    stringsAsFactors = FALSE
  )
}

make_harvest_only <- function(sample_day_tab) {
  key <- paste(sample_day_tab$harvest, sample_day_tab$cohort, sample_day_tab$dose, sep = "__")
  max_day <- ave(sample_day_tab$day, key, FUN = max)
  sample_day_tab[sample_day_tab$day == max_day, , drop = FALSE]
}

write_tsv <- function(tab, path) {
  utils::write.table(tab, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(path)
}

print_signif_df <- function(tab, digits = 6) {
  out <- tab
  is_num <- vapply(out, is.numeric, logical(1))
  for (nm in names(out)[is_num]) {
    out[[nm]] <- signif(out[[nm]], digits = digits)
  }
  print(out)
  invisible(out)
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  resolved <- resolve_seed_dir(argv)
  seed_dir <- resolved$seed_dir
  out_dir <- resolve_path_value(argv$out_dir, getwd()) %||% file.path(seed_dir, "viz", "live_effective_pms")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  params <- read_best_params_map(seed_dir)
  cfg <- read_fit_config(seed_dir)

  ploidy_path <- file.path(seed_dir, "viz", "ploidy_timecourse.tsv")
  o2_path <- file.path(seed_dir, "viz", "predict_burden_vs_o2.tsv")
  viability_path <- file.path(seed_dir, "viz", "functional_curve_ploidy.tsv")
  ploidy_tab <- read_required_tsv(ploidy_path)
  o2_tab <- read_required_tsv(o2_path)
  viability_tab <- read_required_tsv(viability_path)

  needed_ploidy <- c("harvest", "cohort", "dose", "day", "N", "fraction")
  needed_o2 <- c("harvest", "cohort", "dose", "day", "o2_pct")
  needed_viability <- c("N", "viability_after_ms")
  if (!all(needed_ploidy %in% names(ploidy_tab))) {
    stop("ploidy_timecourse.tsv is missing required columns: ", paste(setdiff(needed_ploidy, names(ploidy_tab)), collapse = ", "))
  }
  if (!all(needed_o2 %in% names(o2_tab))) {
    stop("predict_burden_vs_o2.tsv is missing required columns: ", paste(setdiff(needed_o2, names(o2_tab)), collapse = ", "))
  }
  if (!all(needed_viability %in% names(viability_tab))) {
    stop("functional_curve_ploidy.tsv is missing required columns: ", paste(setdiff(needed_viability, names(viability_tab)), collapse = ", "))
  }

  N_UNIT <- as_num(cfg$N_UNIT, 22)
  if (!is.finite(N_UNIT) || N_UNIT <= 0) stop("N_UNIT is missing or invalid in fit_config.rds.")
  N_dip <- 2 * N_UNIT
  ploidy_mode <- as.character(cfg$ploidy_O2_death %||% "ploidy_related")

  merged <- merge(
    ploidy_tab[, needed_ploidy, drop = FALSE],
    o2_tab[, needed_o2, drop = FALSE],
    by = c("harvest", "cohort", "dose", "day"),
    all.x = TRUE
  )
  merged <- merge(
    merged,
    viability_tab[, needed_viability, drop = FALSE],
    by = "N",
    all.x = TRUE
  )
  merged$fraction <- suppressWarnings(as.numeric(merged$fraction))
  merged$N <- suppressWarnings(as.numeric(merged$N))
  merged$o2_pct <- suppressWarnings(as.numeric(merged$o2_pct))
  merged$viability_after_ms <- suppressWarnings(as.numeric(merged$viability_after_ms))
  merged <- merged[
    is.finite(merged$fraction) &
      merged$fraction > 0 &
      is.finite(merged$N) &
      is.finite(merged$o2_pct),
    ,
    drop = FALSE
  ]
  if (!nrow(merged)) {
    stop("No valid live ploidy rows with oxygen values were found after merging ploidy_timecourse.tsv and predict_burden_vs_o2.tsv.")
  }

  p_mis_base <- as_num(params[["p_mis_base"]], 1e-5)
  p_misseg <- as_num(params[["p_misseg"]], 0)
  k_o_mis <- as_num(params[["k_o_mis"]], 1e-12)
  mu_hp <- as_num(params[["mu_hp"]], 0)
  gamma_mu <- as_num(params[["gamma_mu"]], 1)
  O2_crit <- as_num(params[["O2_crit"]], as_num(cfg$O2_crit_init, 1))
  n_O <- as_num(params[["n_O"]], as_num(cfg$n_O_init, 1))

  merged$mu_eff <- compute_mu_eff(
    O2 = merged$o2_pct,
    N = merged$N,
    mu_hp = mu_hp,
    gamma_mu = gamma_mu,
    O2_crit = O2_crit,
    n_O = n_O,
    N_dip = N_dip,
    ploidy_O2_death = ploidy_mode
  )
  merged$p_ms <- compute_p_ms(
    mu_eff = merged$mu_eff,
    p_mis_base = p_mis_base,
    p_misseg = p_misseg,
    k_o_mis = k_o_mis
  )
  merged$viability_after_ms[!is.finite(merged$viability_after_ms)] <- 1.0
  merged$p_ms_retained_proxy <- merged$p_ms * pmax(pmin(merged$viability_after_ms, 1), 0)

  sample_day_tab <- aggregate_sample_day(merged = merged, p_misseg_parameter = p_misseg, N_UNIT = N_UNIT)
  harvest_tab <- make_harvest_only(sample_day_tab)

  overall_tab <- summarise_sample_day_table(sample_day_tab, p_misseg, "all_sample_days")
  harvest_overall_tab <- summarise_sample_day_table(harvest_tab, p_misseg, "harvest_only")
  cohort_all_days_tab <- do.call(
    rbind,
    lapply(
      split(sample_day_tab, sample_day_tab$cohort),
      function(df) summarise_sample_day_table(df, p_misseg, paste0("cohort_", unique(df$cohort), "_all_days"))
    )
  )
  cohort_harvest_tab <- do.call(
    rbind,
    lapply(
      split(harvest_tab, harvest_tab$cohort),
      function(df) summarise_sample_day_table(df, p_misseg, paste0("cohort_", unique(df$cohort), "_harvest_only"))
    )
  )

  context_tab <- data.frame(
    seed_dir = seed_dir,
    seed_id = resolved$seed_id,
    selected_by = resolved$selected_by,
    out_dir = out_dir,
    N_UNIT = N_UNIT,
    N_dip = N_dip,
    ploidy_O2_death = ploidy_mode,
    p_mis_base = p_mis_base,
    p_misseg = p_misseg,
    k_o_mis = k_o_mis,
    mu_hp = mu_hp,
    gamma_mu = gamma_mu,
    O2_crit = O2_crit,
    n_O = n_O,
    stringsAsFactors = FALSE
  )

  write_tsv(context_tab, file.path(out_dir, "live_effective_pms_context.tsv"))
  write_tsv(sample_day_tab, file.path(out_dir, "live_effective_pms_sample_day.tsv"))
  write_tsv(overall_tab, file.path(out_dir, "live_effective_pms_overall.tsv"))
  write_tsv(harvest_overall_tab, file.path(out_dir, "live_effective_pms_harvest_only.tsv"))
  write_tsv(cohort_all_days_tab, file.path(out_dir, "live_effective_pms_cohort_all_days.tsv"))
  write_tsv(cohort_harvest_tab, file.path(out_dir, "live_effective_pms_cohort_harvest_only.tsv"))

  cat("Seed directory:\t", seed_dir, "\n", sep = "")
  cat("Selected by:\t", resolved$selected_by, "\n", sep = "")
  cat("Output dir:\t", out_dir, "\n", sep = "")
  cat("\nOverall live-weighted effective p_ms:\n")
  print_signif_df(overall_tab, 6)
  cat("\nHarvest-only live-weighted effective p_ms:\n")
  print_signif_df(harvest_overall_tab, 6)
  cat("\nCohort-stratified live-weighted effective p_ms (all days):\n")
  print_signif_df(cohort_all_days_tab, 6)
  cat("\nCohort-stratified live-weighted effective p_ms (harvest only):\n")
  print_signif_df(cohort_harvest_tab, 6)
}

if (sys.nframe() == 0) {
  main()
}
