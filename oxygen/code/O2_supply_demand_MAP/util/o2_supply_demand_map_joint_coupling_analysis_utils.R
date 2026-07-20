#!/usr/bin/env Rscript

# Shared, side-effect-free helpers for the joint soft-coupling analyses.
# Executable analysis, visualization, reporting, and HPC orchestration remain in
# their respective layers.

o2jca_null_coalesce <- function(x, y) {
  if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x))) y else x
}

`%||%` <- o2jca_null_coalesce

o2jca_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    bits <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", bits[[1L]], fixed = TRUE)
    out[[key]] <- if (length(bits) > 1L) paste(bits[-1L], collapse = "=") else "TRUE"
  }
  out
}

o2jca_as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1L]])) return(default)
  value <- tolower(trimws(as.character(x[[1L]])))
  if (value %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n", "off")) return(FALSE)
  default
}

o2jca_as_int <- function(x, default = NA_integer_) {
  value <- suppressWarnings(as.integer(x %||% default))
  if (length(value) && is.finite(value[[1L]])) value[[1L]] else default
}

o2jca_as_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x %||% default))
  if (length(value) && is.finite(value[[1L]])) value[[1L]] else default
}

o2jca_split_csv <- function(x, default = character()) {
  if (is.null(x) || !length(x) || is.na(x[[1L]]) || !nzchar(trimws(x[[1L]]))) return(default)
  value <- trimws(strsplit(as.character(x[[1L]]), ",", fixed = TRUE)[[1L]])
  value[nzchar(value)]
}

o2jca_repo_root <- function(workflow_root) {
  normalizePath(file.path(workflow_root, "..", "..", ".."), mustWork = TRUE)
}

o2jca_normalize_path <- function(path, mustWork = FALSE) {
  normalizePath(path, winslash = "/", mustWork = mustWork)
}

o2jca_path_is_within <- function(path, parent) {
  path <- o2jca_normalize_path(path, mustWork = FALSE)
  parent <- sub("/+$", "", o2jca_normalize_path(parent, mustWork = FALSE))
  path == parent | startsWith(path, paste0(parent, "/"))
}

o2jca_default_output_root <- function(result_root, workflow_root) {
  file.path(
    o2jca_repo_root(workflow_root), "oxygen", "results", "analysis",
    "joint_coupling", basename(sub("/+$", "", result_root))
  )
}

o2jca_assert_separate_output_root <- function(result_root, output_root) {
  if (o2jca_path_is_within(output_root, result_root)) {
    stop(
      "output_root must be outside the read-only fitting result_root: ",
      o2jca_normalize_path(output_root, mustWork = FALSE),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

o2jca_pair_short_label <- function(pair_id) {
  meta <- o2jca_pair_metadata(pair_id)
  if (!is.finite(meta$invivo_seed) || !is.finite(meta$invivo_cluster) || !is.finite(meta$invivo_subcluster)) {
    return(pair_id)
  }
  sprintf("C%02dSc%02d / vi%d", meta$invivo_cluster, meta$invivo_subcluster, meta$invivo_seed)
}

o2jca_read_tsv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing required input: ", path, call. = FALSE)
    return(NULL)
  }
  utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

o2jca_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  normalizePath(path, mustWork = TRUE)
}

o2jca_write_manifest <- function(paths, path, artifact_type = "analysis_table") {
  paths <- unique(normalizePath(paths[file.exists(paths)], mustWork = TRUE))
  info <- file.info(paths)
  out <- data.frame(
    artifact_type = artifact_type,
    file = basename(paths),
    path = paths,
    bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  o2jca_write_tsv(out, path)
}

o2jca_seed_number <- function(x) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(x))))
}

o2jca_normalize_seed <- function(x) {
  x <- trimws(as.character(x))
  ifelse(grepl("^[0-9]+$", x), paste0("seed", x), tolower(x))
}

o2jca_parameter_levels <- function() {
  c(
    "O2_crit", "n_O", "alpha_o2", "gamma_growth", "lam_max", "mu_hp",
    "gamma_mu", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
    "buffer_smax", "buffer_beta", "buffer_n_exp"
  )
}

o2jca_process_map <- function() {
  data.frame(
    parameter = o2jca_parameter_levels(),
    primary_process = c(
      "Hypoxia sensing/resource stress", "Hypoxia sensing/resource stress",
      "Growth", "Growth", "Growth", "Hypoxia-associated death",
      "Hypoxia-associated death", "CIN/missegregation", "CIN/missegregation",
      "CIN/missegregation", "Whole-genome doubling",
      "Post-missegregation survival/buffering",
      "Post-missegregation survival/buffering",
      "Post-missegregation survival/buffering"
    ),
    downstream_processes = c(
      "Resource-stress activation; growth and death response",
      "Resource-stress curve steepness; growth and death response",
      "High-chromosome growth damping", "Chromosome dependence of growth penalty",
      "Maximum division rate", "Stress-associated death scale",
      "Chromosome dependence of stress-associated death",
      "Baseline chromosome missegregation", "Stress-induced chromosome missegregation",
      "Stress scale controlling induced missegregation", "WGD generation",
      "Maximum post-missegregation survival", "Ploidy-dependent viability loss",
      "Ploidy exponent of post-missegregation survival"
    ),
    effect_direction = c(
      "Higher threshold activates stress at higher O2", "Higher value sharpens stress transition",
      "Higher value strengthens growth damping", "Higher value strengthens chromosome-dependent growth penalty",
      "Higher value increases maximal division", "Higher value increases stress-associated death",
      "Higher value strengthens chromosome-dependent death penalty",
      "Higher value increases baseline CIN", "Higher value increases stress-induced CIN",
      "Higher value shifts the stress half-saturation scale",
      "Higher value increases WGD events", "Higher value increases maximal survival",
      "Higher value strengthens ploidy-dependent viability loss",
      "Higher value steepens ploidy dependence"
    ),
    code_evidence = c(
      rep("oxygen/data/O2_supply_demand/parameter_table_O2.csv", 14L)
    ),
    mapping_status = "audited_from_parameter_table",
    stringsAsFactors = FALSE
  )
}

o2jca_pair_metadata <- function(pair_id) {
  pattern <- "^fit_joint_(.+)_vi_seed([0-9]+)_C([0-9]+)Sc([0-9]+)_vt_seed([0-9]+)$"
  hit <- regexec(pattern, pair_id)
  values <- regmatches(pair_id, hit)[[1L]]
  if (length(values) != 6L) {
    return(data.frame(
      pair_id = pair_id, warmup_method = NA_character_, invivo_seed = NA_integer_,
      invivo_cluster = NA_integer_, invivo_subcluster = NA_integer_,
      invitro_seed = NA_integer_, stringsAsFactors = FALSE
    ))
  }
  data.frame(
    pair_id = pair_id,
    warmup_method = values[[2L]],
    invivo_seed = as.integer(values[[3L]]),
    invivo_cluster = as.integer(values[[4L]]),
    invivo_subcluster = as.integer(values[[5L]]),
    invitro_seed = as.integer(values[[6L]]),
    stringsAsFactors = FALSE
  )
}

o2jca_discover_pair_dirs <- function(result_root, pair_pattern = "^fit_joint_.*_vt_seed[0-9]+$", max_pairs = NA_integer_) {
  if (!dir.exists(result_root)) stop("Result root does not exist: ", result_root, call. = FALSE)
  candidates <- list.dirs(result_root, recursive = FALSE, full.names = TRUE)
  candidates <- sort(candidates[grepl(pair_pattern, basename(candidates), perl = TRUE)])
  required <- file.path(candidates, "extra_results", "joint_soft_coupling_all.tsv")
  candidates <- candidates[file.exists(required)]
  if (is.finite(max_pairs)) candidates <- head(candidates, max_pairs)
  if (!length(candidates)) stop("No completed joint-pair directories found under: ", result_root, call. = FALSE)
  candidates
}

o2jca_ratio_class <- function(ratio, threshold = 1.1) {
  ratio <- suppressWarnings(as.numeric(ratio))
  lower <- 1 / threshold
  out <- rep("Invalid", length(ratio))
  valid <- is.finite(ratio) & ratio > 0
  out[valid & ratio < lower] <- "ClassA"
  out[valid & ratio >= lower & ratio <= threshold] <- "ClassB"
  out[valid & ratio > threshold] <- "ClassC"
  factor(out, levels = c("ClassA", "ClassB", "ClassC", "Invalid"))
}

o2jca_normalized_entropy <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[is.finite(counts) & counts > 0]
  if (!length(counts) || sum(counts) <= 0) return(NA_real_)
  p <- counts / sum(counts)
  if (length(p) == 1L) return(0)
  -sum(p * log(p)) / log(3)
}

o2jca_safe_stat <- function(x, fun, default = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(default)
  fun(x)
}

o2jca_quantile <- function(x, probability) {
  o2jca_safe_stat(x, function(z) unname(stats::quantile(z, probability, names = FALSE, type = 8)))
}

o2jca_dominant_class <- function(classes) {
  counts <- table(factor(as.character(classes), levels = c("ClassA", "ClassB", "ClassC")))
  if (!sum(counts)) return(c(class = NA_character_, fraction = NA_real_, margin = NA_real_))
  ordered <- sort(as.numeric(counts), decreasing = TRUE)
  winner <- names(counts)[which.max(counts)]
  c(
    class = winner,
    fraction = as.character(ordered[[1L]] / sum(counts)),
    margin = as.character((ordered[[1L]] - ordered[[2L]]) / sum(counts))
  )
}

o2jca_group_split <- function(data, columns) {
  grouping <- do.call(
    interaction,
    c(unclass(data[, columns, drop = FALSE]), list(drop = TRUE, lex.order = TRUE))
  )
  split(data, grouping, drop = TRUE)
}

o2jca_summarize_parameter_group <- function(data, threshold = 1.1) {
  valid <- data[as.character(data$ratio_class) %in% c("ClassA", "ClassB", "ClassC"), , drop = FALSE]
  classes <- factor(as.character(valid$ratio_class), levels = c("ClassA", "ClassB", "ClassC"))
  counts <- table(classes)
  dominant <- o2jca_dominant_class(classes)
  ratio <- suppressWarnings(as.numeric(valid$ratio_vivo_to_vitro))
  log2_ratio <- suppressWarnings(as.numeric(valid$log2_ratio_vivo_to_vitro))
  n_valid <- length(ratio)
  data.frame(
    n_seed_rows = nrow(data),
    n_valid = n_valid,
    n_invalid = nrow(data) - n_valid,
    n_ClassA = unname(counts[["ClassA"]]),
    n_ClassB = unname(counts[["ClassB"]]),
    n_ClassC = unname(counts[["ClassC"]]),
    prop_ClassA = if (n_valid) unname(counts[["ClassA"]]) / n_valid else NA_real_,
    prop_ClassB = if (n_valid) unname(counts[["ClassB"]]) / n_valid else NA_real_,
    prop_ClassC = if (n_valid) unname(counts[["ClassC"]]) / n_valid else NA_real_,
    dominant_class = dominant[["class"]],
    dominant_fraction = as.numeric(dominant[["fraction"]]),
    consensus_margin = as.numeric(dominant[["margin"]]),
    normalized_entropy = o2jca_normalized_entropy(counts),
    union_ClassA = unname(counts[["ClassA"]]) > 0,
    union_ClassB = unname(counts[["ClassB"]]) > 0,
    union_ClassC = unname(counts[["ClassC"]]) > 0,
    intersection_ClassA = n_valid > 0 && unname(counts[["ClassA"]]) == n_valid,
    intersection_ClassB = n_valid > 0 && unname(counts[["ClassB"]]) == n_valid,
    intersection_ClassC = n_valid > 0 && unname(counts[["ClassC"]]) == n_valid,
    stable80_ClassA = n_valid > 0 && unname(counts[["ClassA"]]) / n_valid >= 0.80,
    stable80_ClassB = n_valid > 0 && unname(counts[["ClassB"]]) / n_valid >= 0.80,
    stable80_ClassC = n_valid > 0 && unname(counts[["ClassC"]]) / n_valid >= 0.80,
    stable90_ClassA = n_valid > 0 && unname(counts[["ClassA"]]) / n_valid >= 0.90,
    stable90_ClassB = n_valid > 0 && unname(counts[["ClassB"]]) / n_valid >= 0.90,
    stable90_ClassC = n_valid > 0 && unname(counts[["ClassC"]]) / n_valid >= 0.90,
    stable95_ClassA = n_valid > 0 && unname(counts[["ClassA"]]) / n_valid >= 0.95,
    stable95_ClassB = n_valid > 0 && unname(counts[["ClassB"]]) / n_valid >= 0.95,
    stable95_ClassC = n_valid > 0 && unname(counts[["ClassC"]]) / n_valid >= 0.95,
    ratio_mean = o2jca_safe_stat(ratio, mean),
    ratio_variance = o2jca_safe_stat(ratio, function(z) if (length(z) > 1L) stats::var(z) else 0),
    ratio_median = o2jca_safe_stat(ratio, stats::median),
    ratio_q05 = o2jca_quantile(ratio, 0.05),
    ratio_q25 = o2jca_quantile(ratio, 0.25),
    ratio_q75 = o2jca_quantile(ratio, 0.75),
    ratio_q95 = o2jca_quantile(ratio, 0.95),
    log2_ratio_mean = o2jca_safe_stat(log2_ratio, mean),
    log2_ratio_variance = o2jca_safe_stat(log2_ratio, function(z) if (length(z) > 1L) stats::var(z) else 0),
    log2_ratio_sd = o2jca_safe_stat(log2_ratio, function(z) if (length(z) > 1L) stats::sd(z) else 0),
    log2_ratio_median = o2jca_safe_stat(log2_ratio, stats::median),
    log2_ratio_mad = o2jca_safe_stat(log2_ratio, stats::mad),
    log2_ratio_iqr = o2jca_safe_stat(log2_ratio, stats::IQR),
    geometric_mean_ratio = if (n_valid) 2^mean(log2_ratio) else NA_real_,
    prop_ratio_lt1 = if (n_valid) mean(ratio < 1) else NA_real_,
    prop_ratio_eq1 = if (n_valid) mean(abs(log2_ratio) <= 1e-12) else NA_real_,
    prop_ratio_gt1 = if (n_valid) mean(ratio > 1) else NA_real_,
    mean_abs_log2_distance_from_1 = o2jca_safe_stat(abs(log2_ratio), mean),
    mean_ClassA_threshold_exceedance = o2jca_safe_stat(pmax(0, log2(1 / threshold) - log2_ratio), mean),
    mean_ClassC_threshold_exceedance = o2jca_safe_stat(pmax(0, log2_ratio - log2(threshold)), mean),
    stringsAsFactors = FALSE
  )
}

o2jca_class_long <- function(summary_table, id_columns) {
  rows <- list()
  for (i in seq_len(nrow(summary_table))) {
    for (class_name in c("ClassA", "ClassB", "ClassC")) {
      rows[[length(rows) + 1L]] <- data.frame(
        summary_table[i, id_columns, drop = FALSE],
        ratio_class = class_name,
        n = summary_table[[paste0("n_", class_name)]][[i]],
        proportion = summary_table[[paste0("prop_", class_name)]][[i]],
        union_any_seed = summary_table[[paste0("union_", class_name)]][[i]],
        strict_intersection_all_seeds = summary_table[[paste0("intersection_", class_name)]][[i]],
        stable80 = summary_table[[paste0("stable80_", class_name)]][[i]],
        stable90 = summary_table[[paste0("stable90_", class_name)]][[i]],
        stable95 = summary_table[[paste0("stable95_", class_name)]][[i]],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

o2jca_assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

o2jca_assert_unique_key <- function(data, columns, label) {
  key <- do.call(paste, c(data[, columns, drop = FALSE], sep = "\r"))
  if (anyDuplicated(key)) {
    stop(label, " has duplicate key rows for: ", paste(columns, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

o2jca_objective_strata <- function(seed_summary) {
  seed_summary$objective <- suppressWarnings(as.numeric(seed_summary$objective))
  seed_summary <- seed_summary[is.finite(seed_summary$objective), , drop = FALSE]
  if (!nrow(seed_summary)) return(data.frame())
  seed_summary <- seed_summary[order(seed_summary$objective, o2jca_seed_number(seed_summary$seed)), , drop = FALSE]
  seed_summary$objective_rank_recomputed <- seq_len(nrow(seed_summary))
  seed_summary$objective_delta <- seed_summary$objective - min(seed_summary$objective)
  seed_summary$objective_percentile <- seed_summary$objective_rank_recomputed / nrow(seed_summary)
  rules <- list(
    all_valid = rep(TRUE, nrow(seed_summary)),
    delta_le_2 = seed_summary$objective_delta <= 2,
    delta_le_5 = seed_summary$objective_delta <= 5,
    delta_le_10 = seed_summary$objective_delta <= 10,
    top_10pct = seed_summary$objective_percentile <= 0.10,
    top_25pct = seed_summary$objective_percentile <= 0.25,
    top_50pct = seed_summary$objective_percentile <= 0.50,
    best_seed = seed_summary$objective_rank_recomputed == 1L
  )
  rows <- lapply(names(rules), function(name) {
    out <- seed_summary[rules[[name]], , drop = FALSE]
    if (!nrow(out)) return(NULL)
    out$objective_stratum <- rep(name, nrow(out))
    out
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

o2jca_read_joint_master <- function(path) {
  data <- o2jca_read_tsv(path)
  o2jca_assert_columns(
    data,
    c("pair_id", "seed", "parameter", "ratio_vivo_to_vitro", "log2_ratio_vivo_to_vitro", "ratio_class"),
    "joint soft-coupling master table"
  )
  data
}

o2jca_cramers_v <- function(tab) {
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  n <- sum(tab)
  if (!n || min(dim(tab)) < 2L) return(NA_real_)
  test <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  sqrt(unname(test$statistic) / (n * min(nrow(tab) - 1L, ncol(tab) - 1L)))
}

o2jca_normalized_mutual_information <- function(tab) {
  n <- sum(tab)
  if (!n) return(NA_real_)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  terms <- matrix(0, nrow(tab), ncol(tab))
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (pxy[i, j] > 0 && px[i] > 0 && py[j] > 0) {
        terms[i, j] <- pxy[i, j] * log(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  mi <- sum(terms)
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx <= 0 || hy <= 0) return(0)
  mi / sqrt(hx * hy)
}

o2jca_js_divergence <- function(p, q) {
  p <- as.numeric(p); q <- as.numeric(q)
  if (sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  p <- p / sum(p); q <- q / sum(q); m <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0))
  (kl(p, m) + kl(q, m)) / 2
}

o2jca_bootstrap_pair_mean <- function(values, n_boot = 2000L, seed = 5826L) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (!length(values)) return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  set.seed(seed)
  boot <- replicate(n_boot, mean(sample(values, length(values), replace = TRUE)))
  c(
    mean = mean(values),
    lower = unname(stats::quantile(boot, 0.025, names = FALSE, type = 8)),
    upper = unname(stats::quantile(boot, 0.975, names = FALSE, type = 8))
  )
}
