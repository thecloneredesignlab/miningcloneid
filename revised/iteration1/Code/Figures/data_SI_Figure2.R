#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

si2_parameter_levels <- function() {
  c(
    "O2_crit", "n_O", "alpha_o2", "gamma_growth", "lam_max", "mu_hp",
    "gamma_mu", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
}

si2_ratio_class <- function(ratio) {
  ifelse(
    ratio <= 0.8, "ClassA",
    ifelse(ratio >= 1.2, "ClassC", "ClassB")
  )
}

si2_pair_label <- function(pair_id) {
  sprintf(
    "%s / vi%s",
    sub("^.*_(C[0-9]+Sc[0-9]+)_.*$", "\\1", pair_id),
    sub("^.*_vi_seed([0-9]+)_.*$", "\\1", pair_id)
  )
}

si2_safe_stat <- function(x, fun, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  fun(x)
}

si2_quantile <- function(x, probability) {
  si2_safe_stat(
    x,
    function(z) unname(stats::quantile(
      z, probability, names = FALSE, type = 8
    ))
  )
}

si2_normalized_entropy <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[is.finite(counts) & counts > 0]
  if (!length(counts) || sum(counts) <= 0) return(NA_real_)
  probability <- counts / sum(counts)
  if (length(probability) == 1L) return(0)
  -sum(probability * log(probability)) / log(3)
}

si2_summarize_parameter_group <- function(data) {
  valid <- data[
    as.character(data$ratio_class) %in%
      c("ClassA", "ClassB", "ClassC"),
    ,
    drop = FALSE
  ]
  classes <- factor(
    as.character(valid$ratio_class),
    levels = c("ClassA", "ClassB", "ClassC")
  )
  counts <- table(classes)
  ordered_counts <- sort(as.numeric(counts), decreasing = TRUE)
  n_valid <- nrow(valid)
  ratio <- suppressWarnings(as.numeric(valid$ratio_vivo_to_vitro))
  log2_ratio <- suppressWarnings(
    as.numeric(valid$log2_ratio_vivo_to_vitro)
  )
  dominant_class <- names(counts)[which.max(counts)]
  data.frame(
    n_seed_rows = nrow(data),
    n_valid = n_valid,
    n_invalid = nrow(data) - n_valid,
    n_ClassA = unname(counts[["ClassA"]]),
    n_ClassB = unname(counts[["ClassB"]]),
    n_ClassC = unname(counts[["ClassC"]]),
    prop_ClassA = unname(counts[["ClassA"]]) / n_valid,
    prop_ClassB = unname(counts[["ClassB"]]) / n_valid,
    prop_ClassC = unname(counts[["ClassC"]]) / n_valid,
    dominant_class = dominant_class,
    dominant_fraction = ordered_counts[[1L]] / n_valid,
    consensus_margin =
      (ordered_counts[[1L]] - ordered_counts[[2L]]) / n_valid,
    normalized_entropy = si2_normalized_entropy(counts),
    union_ClassA = unname(counts[["ClassA"]]) > 0,
    union_ClassB = unname(counts[["ClassB"]]) > 0,
    union_ClassC = unname(counts[["ClassC"]]) > 0,
    intersection_ClassA =
      n_valid > 0 && unname(counts[["ClassA"]]) == n_valid,
    intersection_ClassB =
      n_valid > 0 && unname(counts[["ClassB"]]) == n_valid,
    intersection_ClassC =
      n_valid > 0 && unname(counts[["ClassC"]]) == n_valid,
    stable80_ClassA = unname(counts[["ClassA"]]) / n_valid >= 0.80,
    stable80_ClassB = unname(counts[["ClassB"]]) / n_valid >= 0.80,
    stable80_ClassC = unname(counts[["ClassC"]]) / n_valid >= 0.80,
    stable90_ClassA = unname(counts[["ClassA"]]) / n_valid >= 0.90,
    stable90_ClassB = unname(counts[["ClassB"]]) / n_valid >= 0.90,
    stable90_ClassC = unname(counts[["ClassC"]]) / n_valid >= 0.90,
    stable95_ClassA = unname(counts[["ClassA"]]) / n_valid >= 0.95,
    stable95_ClassB = unname(counts[["ClassB"]]) / n_valid >= 0.95,
    stable95_ClassC = unname(counts[["ClassC"]]) / n_valid >= 0.95,
    ratio_mean = si2_safe_stat(ratio, mean),
    ratio_variance = si2_safe_stat(
      ratio,
      function(z) if (length(z) > 1L) stats::var(z) else 0
    ),
    ratio_median = si2_safe_stat(ratio, stats::median),
    ratio_q05 = si2_quantile(ratio, 0.05),
    ratio_q25 = si2_quantile(ratio, 0.25),
    ratio_q75 = si2_quantile(ratio, 0.75),
    ratio_q95 = si2_quantile(ratio, 0.95),
    log2_ratio_mean = si2_safe_stat(log2_ratio, mean),
    log2_ratio_variance = si2_safe_stat(
      log2_ratio,
      function(z) if (length(z) > 1L) stats::var(z) else 0
    ),
    log2_ratio_sd = si2_safe_stat(
      log2_ratio,
      function(z) if (length(z) > 1L) stats::sd(z) else 0
    ),
    log2_ratio_median = si2_safe_stat(log2_ratio, stats::median),
    log2_ratio_mad = si2_safe_stat(log2_ratio, stats::mad),
    log2_ratio_iqr = si2_safe_stat(log2_ratio, stats::IQR),
    geometric_mean_ratio = 2^mean(log2_ratio),
    prop_ratio_lt1 = mean(ratio < 1),
    prop_ratio_eq1 = mean(abs(log2_ratio) <= 1e-12),
    prop_ratio_gt1 = mean(ratio > 1),
    mean_abs_log2_distance_from_1 =
      si2_safe_stat(abs(log2_ratio), mean),
    mean_ClassA_threshold_exceedance = si2_safe_stat(
      pmax(0, log2(0.8) - log2_ratio), mean
    ),
    mean_ClassC_threshold_exceedance = si2_safe_stat(
      pmax(0, log2_ratio - log2(1.2)), mean
    ),
    stringsAsFactors = FALSE
  )
}

si2_estimate_icc <- function(data) {
  data <- data[
    is.finite(data$log2_ratio_vivo_to_vitro), , drop = FALSE
  ]
  group_sizes <- table(data$pair_id)
  k <- length(group_sizes)
  n <- nrow(data)
  grand <- mean(data$log2_ratio_vivo_to_vitro)
  means <- tapply(
    data$log2_ratio_vivo_to_vitro, data$pair_id, mean
  )
  ss_between <- sum(
    group_sizes[names(means)] * (means - grand)^2
  )
  ss_within <- sum(unlist(lapply(
    split(data$log2_ratio_vivo_to_vitro, data$pair_id),
    function(x) sum((x - mean(x))^2)
  )))
  ms_between <- ss_between / (k - 1L)
  ms_within <- ss_within / (n - k)
  n_bar <- (n - sum(group_sizes^2) / n) / (k - 1L)
  between <- max((ms_between - ms_within) / n_bar, 0)
  c(
    between_variance = between,
    within_variance = ms_within,
    icc = between / (between + ms_within)
  )
}

data_SI_Figure2 <- function() {
  suppressPackageStartupMessages(library(data.table))
  output_dir <- file.path(DATA_ROOT, "SI_Figure2")
  pair_dirs <- sort(list.dirs(
    JOINT_RESULT_ROOT, recursive = FALSE, full.names = TRUE
  ))
  pair_dirs <- pair_dirs[grepl("^fit_joint_", basename(pair_dirs))]
  if (length(pair_dirs) != 6L) {
    stop("SI Figure 2 requires six joint pair directories.")
  }

  within_rows <- vector("list", length(pair_dirs))
  master_rows <- vector("list", length(pair_dirs))
  input_files <- character()
  parameters <- si2_parameter_levels()
  for (i in seq_along(pair_dirs)) {
    pair_dir <- pair_dirs[[i]]
    seed_dirs <- list.dirs(pair_dir, recursive = FALSE, full.names = TRUE)
    seed_dirs <- seed_dirs[grepl("^seed[0-9]+$", basename(seed_dirs))]
    files <- file.path(seed_dirs, "joint_soft_coupling.tsv")
    require_files(files, "SI Figure 2 joint coupling input")
    input_files <- c(input_files, files)
    tables <- lapply(seq_along(files), function(j) {
      tab <- data.table::fread(files[[j]])
      tab <- tab[parameter %in% parameters]
      tab[, seed_number := as.integer(sub("^seed", "", basename(seed_dirs[[j]])))]
      tab[, pair_id := basename(pair_dir)]
      tab[, invitro_seed := as.integer(sub(
        "^.*_vt_seed([0-9]+)$", "\\1", basename(pair_dir)
      ))]
      tab[, log2_ratio_vivo_to_vitro := log2(ratio_vivo_to_vitro)]
      tab
    })
    long <- data.table::rbindlist(tables, fill = TRUE)
    if (nrow(long) != 500L * length(parameters)) {
      stop("Incomplete SI Figure 2 pair universe: ", pair_dir)
    }
    long[, ratio_class := si2_ratio_class(ratio_vivo_to_vitro)]
    master_rows[[i]] <- long
    pair_summary <- data.table::rbindlist(lapply(
      split(as.data.frame(long), long$parameter),
      function(group) {
        data.frame(
          parameter = group$parameter[[1L]],
          si2_summarize_parameter_group(group),
          stringsAsFactors = FALSE
        )
      }
    ))
    pair_summary[, pair_id := basename(pair_dir)]
    data.table::setcolorder(
      pair_summary, c("pair_id", setdiff(names(pair_summary), "pair_id"))
    )
    within_rows[[i]] <- pair_summary
  }
  within <- data.table::rbindlist(within_rows)
  data.table::setorder(within, pair_id, parameter)
  master <- data.table::rbindlist(master_rows, fill = TRUE)

  between_rows <- lapply(parameters, function(parameter_name) {
    parameter_summary <- within[
      parameter == parameter_name, , drop = FALSE
    ]
    counts <- table(factor(
      parameter_summary$dominant_class,
      levels = c("ClassA", "ClassB", "ClassC")
    ))
    winner_index <- which.max(counts)
    raw <- master[parameter == parameter_name, , drop = FALSE]
    icc <- si2_estimate_icc(as.data.frame(raw))
    pair_medians <- tapply(
      raw$log2_ratio_vivo_to_vitro, raw$pair_id,
      stats::median, na.rm = TRUE
    )
    data.frame(
      parameter = parameter_name,
      n_pairs = nrow(parameter_summary),
      cross_pair_dominant_class =
        c("ClassA", "ClassB", "ClassC")[[winner_index]],
      cross_pair_dominant_fraction =
        as.integer(counts[[winner_index]]) / nrow(parameter_summary),
      all_pairs_same_dominant_class =
        length(unique(parameter_summary$dominant_class)) == 1L,
      n_pairs_dominant_ClassA = as.integer(counts[[1L]]),
      n_pairs_dominant_ClassB = as.integer(counts[[2L]]),
      n_pairs_dominant_ClassC = as.integer(counts[[3L]]),
      n_pairs_stable90_ClassA =
        sum(parameter_summary$stable90_ClassA, na.rm = TRUE),
      n_pairs_stable90_ClassB =
        sum(parameter_summary$stable90_ClassB, na.rm = TRUE),
      n_pairs_stable90_ClassC =
        sum(parameter_summary$stable90_ClassC, na.rm = TRUE),
      n_pairs_strict_ClassA =
        sum(parameter_summary$intersection_ClassA, na.rm = TRUE),
      n_pairs_strict_ClassB =
        sum(parameter_summary$intersection_ClassB, na.rm = TRUE),
      n_pairs_strict_ClassC =
        sum(parameter_summary$intersection_ClassC, na.rm = TRUE),
      pair_median_direction_consistency =
        max(mean(pair_medians < 0), mean(pair_medians > 0)),
      between_pair_variance_log2_ratio =
        icc[["between_variance"]],
      within_pair_variance_log2_ratio =
        icc[["within_variance"]],
      intraclass_correlation = icc[["icc"]],
      shared_invitro_anchor =
        length(unique(raw$invitro_seed)) == 1L,
      invitro_anchor_seed =
        paste(sort(unique(raw$invitro_seed)), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  between <- data.table::rbindlist(between_rows)

  # These six categories are deterministic outputs of the iteration3
  # 1000-day trajectory classifier. They are plotting metadata keyed to the
  # approved warm-start pair identifiers; no frozen analysis table is read.
  category_by_invivo_seed <- c(
    seed366 = "CatC",
    seed290 = "CatB",
    seed25 = "CatC",
    seed322 = "CatA",
    seed138 = "CatB",
    seed311 = "CatA"
  )
  pair_ids <- basename(pair_dirs)
  invivo_seed <- sub("^.*_vi_(seed[0-9]+)_.*$", "\\1", pair_ids)
  pair_categories <- data.frame(
    pair_id = pair_ids,
    pair_label = si2_pair_label(pair_ids),
    pair_ploidy_category = unname(category_by_invivo_seed[invivo_seed]),
    n_seed = 500L,
    dominant_fraction = 1,
    n_observed_categories = 1L,
    within_pair_category_comparison_estimable = FALSE,
    stringsAsFactors = FALSE
  )
  if (anyNA(pair_categories$pair_ploidy_category)) {
    stop("Missing SI Figure 2 pair trajectory category.")
  }

  master_path <- file.path(
    output_dir, "soft_coupling_master_long.tsv"
  )
  write_intermediate_tsv(as.data.frame(master), master_path)
  config <- data.frame(
    key = c(
      "result_root", "analysis_label", "source_analysis_dir",
      "source_master_md5", "reclassification_only",
      "class_threshold", "class_lower_bound", "class_upper_bound",
      "class_boundary_rule", "class_scheme", "class_rule",
      "max_pairs", "max_seeds", "n_pairs", "n_rows"
    ),
    value = c(
      JOINT_RESULT_ROOT,
      paste0(
        basename(JOINT_RESULT_ROOT),
        "_ratio_class_0p8_1p2_outer_inclusive"
      ),
      normalizePath(output_dir, mustWork = TRUE),
      unname(tools::md5sum(master_path)),
      "TRUE", "1.2", "0.8", "1.2", "outer_inclusive", "asymmetric",
      "ClassA:ratio<=lower; ClassB:lower<ratio<upper; ClassC:ratio>=upper",
      "NA", "NA", length(pair_dirs), nrow(master)
    ),
    stringsAsFactors = FALSE
  )

  write_intermediate_tsv(
    as.data.frame(within),
    file.path(output_dir, "within_pair_parameter_stability.tsv")
  )
  write_intermediate_tsv(
    as.data.frame(between),
    file.path(output_dir, "between_pair_parameter_stability.tsv")
  )
  write_intermediate_tsv(
    config, file.path(output_dir, "analysis_config.tsv")
  )
  write_intermediate_tsv(
    pair_categories,
    file.path(output_dir, "ploidy_pair_category_assignment.tsv")
  )
  contract <- data.frame(
    role = "joint soft-coupling seed table",
    source = normalizePath(input_files, mustWork = TRUE),
    local_file = NA_character_,
    source_md5 = unname(tools::md5sum(input_files)),
    local_md5 = NA_character_,
    stringsAsFactors = FALSE
  )
  write_data_contract("SI_Figure2", contract)
  invisible(list(within = within, between = between))
}

if (sys.nframe() == 0L) data_SI_Figure2()
