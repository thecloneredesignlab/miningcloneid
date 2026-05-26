#!/usr/bin/env Rscript

usage <- function() {
  cat(
"Usage:
  Rscript derive_joint_anchor_from_single_fits.R \\
    --invivo_run_dir=/path/to/invivo_run \\
    --invitro_run_dir=/path/to/invitro_run \\
    --base_config=/path/to/O2G_supply_demand.yaml \\
    --out_dir=/path/to/joint_run_dir \\
    --out_config_yaml=/path/to/joint_anchor_config.yaml \\
    --run_prefix=fit_joint_run_prefix \\
    [--anchor_mode=auto|manual] \\
    [--manual_o2_grid=0,0.1,0.2,0.5,1,2] \\
    [--manual_n_grid=44,66,88] \\
    [--auto_o2_bin_pct=0.1] \\
    [--auto_n_bin=1] \\
    [--auto_n_quantiles=0.1,0.5,0.9] \\
    [--auto_max_dynamic=Inf] \\
    [--auto_include_reference=TRUE] \\
    [--top_n=3] \\
    [--extra_results_script=/path/to/extra_results.R] \\
    [--run_extra_results=TRUE]

Description:
  Selects best top-N seeds from finished in vivo and in vitro single-fit runs,
  writes the joint fixed-anchor table, and creates a joint config snapshot.

Selection rules:
  in vivo: finite objective with pred1000_2N > threshold_N and pred1000_4N > threshold_N,
           ranked by minimum objective.
  in vitro: finite objective, ranked by minimum objective.

Anchor modes:
  auto:   derives fixed joint anchors from top in vivo seed trajectories by
          combining O2 timecourses with chromosome-count distributions.
  manual: uses --manual_o2_grid and --manual_n_grid exactly.

Outputs:
  joint_anchor.tsv
  joint_anchor_summary.md
  selected_single_fit_top3.tsv
  joint_anchor_config.yaml
", sep = "")
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    pos <- regexpr("=", kv, fixed = TRUE)
    if (pos < 0) {
      out[[kv]] <- TRUE
    } else {
      key <- substr(kv, 1, pos - 1)
      val <- substr(kv, pos + 1, nchar(kv))
      out[[key]] <- val
    }
  }
  out
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L || is.na(a[[1]])) b else a
}

as_chr <- function(x, default = "") {
  val <- as.character(x %||% default)
  if (!length(val) || !nzchar(trimws(val[[1]]))) default else val[[1]]
}

as_num <- function(x, default = NA_real_) {
  val <- suppressWarnings(as.numeric(x %||% default))
  if (!length(val) || !is.finite(val[[1]])) default else val[[1]]
}

as_num_allow_inf <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  txt <- tolower(trimws(as.character(x[[1]])))
  if (txt %in% c("inf", "+inf", "infinity", "+infinity")) return(Inf)
  if (txt %in% c("-inf", "-infinity")) return(-Inf)
  val <- suppressWarnings(as.numeric(x[[1]]))
  if (!length(val) || is.na(val[[1]])) default else val[[1]]
}

as_int <- function(x, default = NA_integer_) {
  val <- suppressWarnings(as.integer(x %||% default))
  if (!length(val) || !is.finite(val[[1]])) default else val[[1]]
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  val <- tolower(trimws(as.character(x[[1]])))
  if (val %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (val %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)
  default
}

split_num <- function(x, default) {
  txt <- trimws(as_chr(x, paste(default, collapse = ",")))
  if (!nzchar(txt)) return(as.numeric(default))
  vals <- suppressWarnings(as.numeric(trimws(strsplit(txt, "[,;[:space:]]+")[[1]])))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) as.numeric(default) else vals
}

as_num_vec <- function(x, default) {
  if (is.null(x) || !length(x)) return(as.numeric(default))
  if (is.character(x) && length(x) == 1L) {
    vals <- suppressWarnings(as.numeric(trimws(strsplit(x, "[,;[:space:]]+")[[1]])))
  } else {
    vals <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  }
  vals <- vals[is.finite(vals)]
  if (!length(vals)) as.numeric(default) else vals
}

clip_num <- function(x, lower = -Inf, upper = Inf) {
  pmin(pmax(x, lower), upper)
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }
  normalizePath("oxygen/code/O2G_supply_demand_MAP/analysis/derive_joint_anchor_from_single_fits.R", mustWork = FALSE)
}

run_extra_results_if_needed <- function(run_dir, extra_results_script, force = FALSE) {
  summary_path <- file.path(run_dir, "extra_results", "seed_summary.tsv")
  if (!isTRUE(force) && file.exists(summary_path)) {
    return(summary_path)
  }
  if (!file.exists(extra_results_script)) {
    stop("extra_results_script does not exist: ", extra_results_script, call. = FALSE)
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(run_dir, "extra_results_run.log")
  args <- c(extra_results_script, paste0("--run_dir=", run_dir))
  status <- system2("Rscript", args = args, stdout = log_path, stderr = log_path)
  if (!identical(status, 0L) || !file.exists(summary_path)) {
    detail <- if (file.exists(log_path)) paste(utils::tail(readLines(log_path, warn = FALSE), 40L), collapse = "\n") else "(missing log)"
    stop("extra_results failed for ", run_dir, "\n", detail, call. = FALSE)
  }
  summary_path
}

required_path <- function(x, label, must_work = TRUE) {
  val <- as_chr(x)
  if (!nzchar(val)) stop(label, " is required.", call. = FALSE)
  normalizePath(val, mustWork = must_work)
}

seed_number <- function(seed) {
  suppressWarnings(as.integer(sub("^seed", "", as.character(seed))))
}

normalize_seed_label <- function(seed) {
  txt <- as.character(seed)
  ifelse(grepl("^seed[0-9]+$", txt), txt, paste0("seed", seed_number(txt)))
}

finite_col <- function(tab, candidates) {
  hits <- intersect(candidates, names(tab))
  if (!length(hits)) return(rep(NA_real_, nrow(tab)))
  vals <- suppressWarnings(as.numeric(tab[[hits[[1]]]]))
  vals
}

select_top_seeds <- function(summary_path, scene, top_n = 3L, threshold_N = 44) {
  tab <- read_tsv(summary_path)
  if (!nrow(tab)) stop("summary has no rows: ", summary_path, call. = FALSE)
  if (!"seed" %in% names(tab)) stop("summary missing seed column: ", summary_path, call. = FALSE)
  objective <- finite_col(tab, c("objective", "objective_total", "optimizer_local_objective", "optimizer_deoptim_objective"))
  tab$selection_objective <- objective
  tab$seed <- normalize_seed_label(tab$seed)
  tab$seed_number <- seed_number(tab$seed)

  if (identical(scene, "in_vivo")) {
    pred2 <- finite_col(tab, c("pred1000_2N", "predict1000_2N", "pred_1000_2N"))
    pred4 <- finite_col(tab, c("pred1000_4N", "predict1000_4N", "pred_1000_4N"))
    tab$selection_pred1000_2N <- pred2
    tab$selection_pred1000_4N <- pred4
    tab$eligible <- is.finite(tab$selection_objective) &
      is.finite(pred2) & is.finite(pred4) &
      pred2 > threshold_N & pred4 > threshold_N
    tab$selection_rule <- paste0("objective finite; pred1000_2N>", threshold_N, "; pred1000_4N>", threshold_N)
  } else {
    tab$selection_pred1000_2N <- NA_real_
    tab$selection_pred1000_4N <- NA_real_
    tab$eligible <- is.finite(tab$selection_objective)
    tab$selection_rule <- "objective finite"
  }

  eligible <- tab[tab$eligible %in% TRUE, , drop = FALSE]
  if (!nrow(eligible)) {
    stop("No eligible ", scene, " seed found in ", summary_path, call. = FALSE)
  }
  eligible <- eligible[order(eligible$selection_objective, eligible$seed_number, eligible$seed), , drop = FALSE]
  eligible <- head(eligible, top_n)
  data.frame(
    scene = scene,
    rank = seq_len(nrow(eligible)),
    seed = eligible$seed,
    seed_number = eligible$seed_number,
    objective = eligible$selection_objective,
    pred1000_2N = eligible$selection_pred1000_2N,
    pred1000_4N = eligible$selection_pred1000_4N,
    selection_rule = eligible$selection_rule,
    stringsAsFactors = FALSE
  )
}

anchor_weights <- function(anchor_df, low_o2_weight = 2, zero_o2_priority_weight = 4, n_dip = 44) {
  w <- rep(1, nrow(anchor_df))
  w[is.finite(anchor_df$O2) & anchor_df$O2 <= 0.2] <- low_o2_weight
  w[abs(anchor_df$O2) < 1e-12 & anchor_df$N %in% c(n_dip, 2 * n_dip)] <- zero_o2_priority_weight
  w
}

weighted_quantile <- function(x, w, probs) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  probs <- as.numeric(probs)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / pmax(sum(w), 1e-12)
  vapply(probs, function(p) {
    p <- min(max(as.numeric(p), 0), 1)
    idx <- which(cw >= p)[1]
    if (is.na(idx)) tail(x, 1) else x[[idx]]
  }, numeric(1))
}

first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits)) hits[[1]] else ""
}

bind_rows_fill <- function(...) {
  tabs <- list(...)
  tabs <- tabs[vapply(tabs, is.data.frame, logical(1))]
  tabs <- tabs[vapply(tabs, nrow, integer(1)) > 0]
  if (!length(tabs)) return(data.frame())
  cols <- unique(unlist(lapply(tabs, names), use.names = FALSE))
  tabs <- lapply(tabs, function(tab) {
    missing <- setdiff(cols, names(tab))
    for (col in missing) tab[[col]] <- NA
    tab[, cols, drop = FALSE]
  })
  do.call(rbind, tabs)
}

seed_viz_file <- function(seed_dir, candidates) {
  first_existing(file.path(seed_dir, "viz", candidates))
}

required_cols <- function(tab, cols, label) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) {
    stop(label, " is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

anchor_group_key <- function(seed, harvest, cohort, dose, day) {
  paste(
    seed,
    as.character(harvest),
    as.character(cohort),
    formatC(suppressWarnings(as.numeric(dose)), digits = 12, format = "fg"),
    formatC(suppressWarnings(as.numeric(day)), digits = 12, format = "fg"),
    sep = "\t"
  )
}

anchor_scenario_key <- function(seed, harvest, cohort, dose) {
  paste(
    seed,
    as.character(harvest),
    as.character(cohort),
    formatC(suppressWarnings(as.numeric(dose)), digits = 12, format = "fg"),
    sep = "\t"
  )
}

seed_trajectory_anchor_candidates <- function(invivo_run_dir, seed, quantile_probs) {
  seed_dir <- file.path(invivo_run_dir, seed)
  if (!dir.exists(seed_dir)) {
    stop("Selected in vivo seed directory does not exist: ", seed_dir, call. = FALSE)
  }
  burden_path <- seed_viz_file(
    seed_dir,
    c("predict_burden_0_1000day.tsv", "predict_burden_0_300day.tsv", "predict_burden_0_100day.tsv")
  )
  ploidy_path <- seed_viz_file(
    seed_dir,
    c("predict_ploidy_0_1000day.tsv", "predict_ploidy_0_300day.tsv", "predict_ploidy_0_100day.tsv")
  )
  if (!nzchar(burden_path) || !nzchar(ploidy_path)) {
    stop(
      "Auto anchor requires per-seed viz predictions for ", seed,
      ". Missing ",
      if (!nzchar(burden_path)) "predict_burden_0_*day.tsv" else "",
      if (!nzchar(burden_path) && !nzchar(ploidy_path)) " and " else "",
      if (!nzchar(ploidy_path)) "predict_ploidy_0_*day.tsv" else "",
      " under ", file.path(seed_dir, "viz"),
      call. = FALSE
    )
  }

  burden <- read_tsv(burden_path)
  ploidy <- read_tsv(ploidy_path)
  required_cols(burden, c("harvest", "cohort", "dose", "day", "pred_o2_pct"), burden_path)
  required_cols(ploidy, c("harvest", "cohort", "dose", "day", "N", "fraction"), ploidy_path)

  burden <- burden[, c("harvest", "cohort", "dose", "day", "pred_o2_pct"), drop = FALSE]
  names(burden)[names(burden) == "pred_o2_pct"] <- "O2"
  burden$dose <- suppressWarnings(as.numeric(burden$dose))
  burden$day <- suppressWarnings(as.numeric(burden$day))
  burden$O2 <- suppressWarnings(as.numeric(burden$O2))
  burden <- burden[is.finite(burden$dose) & is.finite(burden$day) & is.finite(burden$O2), , drop = FALSE]

  ploidy <- ploidy[, c("harvest", "cohort", "dose", "day", "N", "fraction"), drop = FALSE]
  ploidy$dose <- suppressWarnings(as.numeric(ploidy$dose))
  ploidy$day <- suppressWarnings(as.numeric(ploidy$day))
  ploidy$N <- suppressWarnings(as.numeric(ploidy$N))
  ploidy$fraction <- suppressWarnings(as.numeric(ploidy$fraction))
  ploidy <- ploidy[
    is.finite(ploidy$dose) & is.finite(ploidy$day) &
      is.finite(ploidy$N) & is.finite(ploidy$fraction),
    ,
    drop = FALSE
  ]
  if (!nrow(burden) || !nrow(ploidy)) return(data.frame())

  merged <- merge(
    ploidy,
    burden,
    by = c("harvest", "cohort", "dose", "day"),
    all = FALSE,
    sort = FALSE
  )
  merged <- merged[is.finite(merged$O2), , drop = FALSE]
  if (!nrow(merged)) return(data.frame())

  keys <- anchor_group_key(seed, merged$harvest, merged$cohort, merged$dose, merged$day)
  split_idx <- split(seq_len(nrow(merged)), keys)
  rows <- vector("list", length(split_idx))
  row_i <- 0L
  q_labels <- paste0("q", formatC(100 * quantile_probs, format = "fg", digits = 4))

  for (idx in split_idx) {
    sub <- merged[idx, , drop = FALSE]
    w <- pmax(as.numeric(sub$fraction), 0)
    if (!(sum(w, na.rm = TRUE) > 0)) next
    n_vals <- as.numeric(sub$N)
    o2 <- as.numeric(sub$O2[[which(is.finite(sub$O2))[1]]])
    if (!is.finite(o2)) next
    mean_N <- sum(n_vals * w, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-12)
    mode_N <- n_vals[which.max(w)]
    q_N <- weighted_quantile(n_vals, w, quantile_probs)
    anchor_N <- c(mean_N, mode_N, q_N)
    anchor_stat <- c("mean", "mode", q_labels)

    row_i <- row_i + 1L
    rows[[row_i]] <- data.frame(
      seed = seed,
      harvest = as.character(sub$harvest[[1]]),
      cohort = as.character(sub$cohort[[1]]),
      dose = as.numeric(sub$dose[[1]]),
      day = as.numeric(sub$day[[1]]),
      O2 = rep(o2, length(anchor_N)),
      N = as.numeric(anchor_N),
      anchor_stat = anchor_stat,
      anchor_source = "top_seed_trajectory",
      anchor_role = "trajectory",
      candidate_weight = 1.0,
      scenario_key = anchor_scenario_key(seed, sub$harvest[[1]], sub$cohort[[1]], sub$dose[[1]]),
      source_burden_file = burden_path,
      source_ploidy_file = ploidy_path,
      stringsAsFactors = FALSE
    )
  }
  if (!row_i) return(data.frame())
  rows <- rows[seq_len(row_i)]
  out <- do.call(rbind, rows)
  if (!is.data.frame(out) || !nrow(out)) return(data.frame())

  priority <- vector("list", 0L)
  for (scenario in unique(out$scenario_key)) {
    sub <- out[out$scenario_key == scenario, , drop = FALSE]
    if (!nrow(sub)) next
    min_o2_day <- sub$day[which.min(sub$O2)]
    endpoint_day <- max(sub$day, na.rm = TRUE)
    for (role in c("min_o2", "endpoint")) {
      day_use <- if (identical(role, "min_o2")) min_o2_day else endpoint_day
      pick <- sub[is.finite(sub$day) & abs(sub$day - day_use) < 1e-8, , drop = FALSE]
      if (!nrow(pick)) next
      pick$anchor_role <- role
      pick$candidate_weight <- 2.0
      priority[[length(priority) + 1L]] <- pick
    }
  }
  if (length(priority)) out <- rbind(out, do.call(rbind, priority))
  out
}

collapse_anchor_candidates <- function(anchor_df, o2_bin, n_bin, n_min, n_max, source_default = "top_seed_trajectory") {
  if (!is.data.frame(anchor_df) || !nrow(anchor_df)) return(data.frame())
  anchor_df$O2 <- suppressWarnings(as.numeric(anchor_df$O2))
  anchor_df$N <- suppressWarnings(as.numeric(anchor_df$N))
  anchor_df <- anchor_df[is.finite(anchor_df$O2) & is.finite(anchor_df$N), , drop = FALSE]
  if (!nrow(anchor_df)) return(anchor_df)
  anchor_df$O2 <- round(clip_num(anchor_df$O2, 0, 100) / o2_bin) * o2_bin
  anchor_df$N <- round(anchor_df$N / n_bin) * n_bin
  anchor_df$N <- clip_num(anchor_df$N, n_min, n_max)
  if (!"candidate_weight" %in% names(anchor_df)) anchor_df$candidate_weight <- 1.0
  anchor_df$candidate_weight <- suppressWarnings(as.numeric(anchor_df$candidate_weight))
  anchor_df$candidate_weight[!is.finite(anchor_df$candidate_weight) | anchor_df$candidate_weight <= 0] <- 1.0
  if (!"anchor_source" %in% names(anchor_df)) anchor_df$anchor_source <- source_default
  if (!"anchor_role" %in% names(anchor_df)) anchor_df$anchor_role <- source_default

  keys <- paste(formatC(anchor_df$O2, digits = 12, format = "fg"), formatC(anchor_df$N, digits = 12, format = "fg"), sep = "\t")
  split_idx <- split(seq_len(nrow(anchor_df)), keys)
  rows <- lapply(split_idx, function(idx) {
    sub <- anchor_df[idx, , drop = FALSE]
    scenario_vals <- if ("scenario_key" %in% names(sub)) {
      vals <- as.character(sub$scenario_key)
      unique(vals[!is.na(vals) & nzchar(vals)])
    } else {
      character()
    }
    data.frame(
      O2 = as.numeric(sub$O2[[1]]),
      N = as.numeric(sub$N[[1]]),
      anchor_source = paste(sort(unique(as.character(sub$anchor_source))), collapse = "+"),
      anchor_role = paste(sort(unique(as.character(sub$anchor_role))), collapse = "+"),
      anchor_weight = sqrt(sum(sub$candidate_weight, na.rm = TRUE)),
      anchor_n_points = nrow(sub),
      anchor_n_scenarios = if (length(scenario_vals)) length(scenario_vals) else NA_integer_,
      anchor_min_day = if ("day" %in% names(sub)) min(suppressWarnings(as.numeric(sub$day)), na.rm = TRUE) else NA_real_,
      anchor_max_day = if ("day" %in% names(sub)) max(suppressWarnings(as.numeric(sub$day)), na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$anchor_weight[!is.finite(out$anchor_weight) | out$anchor_weight <= 0] <- 1.0
  out$anchor_min_day[!is.finite(out$anchor_min_day)] <- NA_real_
  out$anchor_max_day[!is.finite(out$anchor_max_day)] <- NA_real_
  out[order(out$O2, out$N), , drop = FALSE]
}

limit_dynamic_anchors <- function(anchor_df, max_dynamic) {
  if (!is.finite(max_dynamic) || max_dynamic <= 0 || nrow(anchor_df) <= max_dynamic) return(anchor_df)
  ord <- order(-as.numeric(anchor_df$anchor_weight), anchor_df$O2, anchor_df$N)
  keep <- ord[seq_len(as.integer(max_dynamic))]
  out <- anchor_df[keep, , drop = FALSE]
  out[order(out$O2, out$N), , drop = FALSE]
}

plus_union <- function(a, b) {
  vals <- unique(c(strsplit(as.character(a), "\\+", fixed = FALSE)[[1]], as.character(b)))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  paste(sort(vals), collapse = "+")
}

add_reference_anchors <- function(anchors, n_unit) {
  ref <- data.frame(
    O2 = c(0, 0),
    N = c(2 * n_unit, 4 * n_unit),
    anchor_source = "reference_boundary",
    anchor_role = "zero_o2_reference",
    anchor_weight = 1.0,
    anchor_n_points = NA_integer_,
    anchor_n_scenarios = NA_integer_,
    anchor_min_day = NA_real_,
    anchor_max_day = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(ref))) {
    hit <- which(abs(as.numeric(anchors$O2) - ref$O2[[i]]) < 1e-12 &
      abs(as.numeric(anchors$N) - ref$N[[i]]) < 1e-12)
    if (length(hit)) {
      j <- hit[[1]]
      anchors$anchor_source[[j]] <- plus_union(anchors$anchor_source[[j]], ref$anchor_source[[i]])
      anchors$anchor_role[[j]] <- plus_union(anchors$anchor_role[[j]], ref$anchor_role[[i]])
      anchors$anchor_weight[[j]] <- sqrt(as.numeric(anchors$anchor_weight[[j]])^2 + ref$anchor_weight[[i]]^2)
      if ("anchor_n_points" %in% names(anchors) && is.finite(suppressWarnings(as.numeric(anchors$anchor_n_points[[j]])))) {
        anchors$anchor_n_points[[j]] <- as.integer(anchors$anchor_n_points[[j]]) + 1L
      }
    } else {
      anchors <- bind_rows_fill(anchors, ref[i, , drop = FALSE])
    }
  }
  anchors[order(anchors$O2, anchors$N), , drop = FALSE]
}

derive_auto_anchor <- function(invivo_run_dir, selected, cfg, argv) {
  invivo_selected <- selected[selected$scene == "in_vivo", , drop = FALSE]
  if (!nrow(invivo_selected)) stop("No selected in vivo seeds available for auto anchor derivation.", call. = FALSE)

  quantile_probs <- as_num_vec(
    argv$auto_n_quantiles %||% argv$joint_composite_anchor_n_quantiles %||% cfg[["joint_composite_anchor_n_quantiles"]],
    c(0.1, 0.5, 0.9)
  )
  quantile_probs <- sort(unique(clip_num(quantile_probs, 0, 1)))
  if (!length(quantile_probs)) quantile_probs <- c(0.1, 0.5, 0.9)
  o2_bin <- as_num(
    argv$auto_o2_bin_pct %||% argv$joint_composite_anchor_o2_bin_pct %||% cfg[["joint_composite_anchor_o2_bin_pct"]],
    0.1
  )
  n_bin <- as_num(
    argv$auto_n_bin %||% argv$joint_composite_anchor_n_bin %||% cfg[["joint_composite_anchor_n_bin"]],
    1
  )
  if (!is.finite(o2_bin) || o2_bin <= 0) o2_bin <- 0.1
  if (!is.finite(n_bin) || n_bin <= 0) n_bin <- 1
  max_dynamic <- as_num_allow_inf(
    argv$auto_max_dynamic %||% argv$joint_composite_anchor_max_dynamic %||% cfg[["joint_composite_anchor_max_dynamic"]],
    Inf
  )
  include_reference <- as_bool(
    argv$auto_include_reference %||% argv$joint_composite_anchor_include_reference %||% cfg[["joint_composite_anchor_include_reference"]],
    TRUE
  )
  n_unit <- as_num(cfg[["N_UNIT"]], 22)
  if (!is.finite(n_unit) || n_unit <= 0) n_unit <- 22
  n_min <- as_num(cfg[["N_MIN"]], n_unit)
  n_max <- as_num(cfg[["N_MAX"]], 7 * n_unit)
  if (!is.finite(n_min)) n_min <- n_unit
  if (!is.finite(n_max) || n_max < n_min) n_max <- 7 * n_unit

  candidates <- lapply(invivo_selected$seed, function(seed) {
    seed_trajectory_anchor_candidates(invivo_run_dir, seed, quantile_probs = quantile_probs)
  })
  candidates <- candidates[vapply(candidates, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (!length(candidates)) {
    stop("Auto anchor derivation found no trajectory candidates in selected in vivo seeds.", call. = FALSE)
  }
  candidates <- do.call(rbind, candidates)
  anchors <- collapse_anchor_candidates(
    candidates,
    o2_bin = o2_bin,
    n_bin = n_bin,
    n_min = n_min,
    n_max = n_max,
    source_default = "top_seed_trajectory"
  )
  anchors <- limit_dynamic_anchors(anchors, max_dynamic)
  if (isTRUE(include_reference)) {
    anchors <- add_reference_anchors(anchors, n_unit = n_unit)
  }

  attr(anchors, "auto_anchor_settings") <- list(
    o2_bin = o2_bin,
    n_bin = n_bin,
    n_quantiles = quantile_probs,
    max_dynamic = max_dynamic,
    include_reference = include_reference,
    n_candidates = nrow(candidates),
    n_selected_invivo_seeds = nrow(invivo_selected)
  )
  anchors
}

absolutize_if_present <- function(cfg, key, config_dir) {
  val <- cfg[[key]]
  if (is.null(val) || !length(val) || is.na(val[[1]]) || !nzchar(as.character(val[[1]]))) return(cfg)
  txt <- as.character(val[[1]])
  if (!grepl("^/", txt)) txt <- normalizePath(file.path(config_dir, txt), mustWork = FALSE)
  cfg[[key]] <- txt
  cfg
}

main <- function(argv = parse_args(commandArgs(trailingOnly = TRUE))) {
  if (isTRUE(argv$help)) {
    usage()
    return(invisible(0L))
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The R package 'yaml' is required.", call. = FALSE)
  }

  script_dir <- dirname(script_path())
  workflow_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
  oxygen_root <- normalizePath(file.path(workflow_root, "..", ".."), mustWork = FALSE)
  project_root <- normalizePath(file.path(oxygen_root, ".."), mustWork = FALSE)
  default_extra <- file.path(script_dir, "extra_results.R")

  invivo_run_dir <- required_path(argv$invivo_run_dir, "invivo_run_dir")
  invitro_run_dir <- required_path(argv$invitro_run_dir, "invitro_run_dir")
  base_config <- required_path(argv$base_config, "base_config")
  out_dir <- required_path(argv$out_dir, "out_dir", must_work = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_anchor_tsv <- normalizePath(as_chr(argv$out_anchor_tsv, file.path(out_dir, "joint_anchor.tsv")), mustWork = FALSE)
  out_summary_md <- normalizePath(as_chr(argv$out_summary_md, file.path(out_dir, "joint_anchor_summary.md")), mustWork = FALSE)
  out_selected_tsv <- normalizePath(as_chr(argv$out_selected_tsv, file.path(out_dir, "selected_single_fit_top3.tsv")), mustWork = FALSE)
  out_config_yaml <- normalizePath(as_chr(argv$out_config_yaml, file.path(out_dir, "joint_anchor_config.yaml")), mustWork = FALSE)
  extra_results_script <- normalizePath(as_chr(argv$extra_results_script, default_extra), mustWork = FALSE)
  run_extra <- as_bool(argv$run_extra_results, TRUE)
  force_extra <- as_bool(argv$force_extra_results, FALSE)
  top_n <- as_int(argv$top_n, 3L)
  threshold_N <- as_num(argv$threshold_N, 44)
  anchor_mode <- tolower(trimws(as_chr(argv$anchor_mode, "auto")))
  if (!anchor_mode %in% c("auto", "manual")) {
    stop("anchor_mode must be auto or manual, got: ", anchor_mode, call. = FALSE)
  }
  config_dir <- dirname(base_config)
  cfg <- yaml::read_yaml(base_config)

  invivo_summary <- file.path(invivo_run_dir, "extra_results", "seed_summary.tsv")
  invitro_summary <- file.path(invitro_run_dir, "extra_results", "seed_summary.tsv")
  if (isTRUE(run_extra)) {
    invivo_summary <- run_extra_results_if_needed(invivo_run_dir, extra_results_script, force = force_extra)
    invitro_summary <- run_extra_results_if_needed(invitro_run_dir, extra_results_script, force = force_extra)
  } else {
    if (!file.exists(invivo_summary)) stop("Missing in vivo extra_results seed_summary: ", invivo_summary, call. = FALSE)
    if (!file.exists(invitro_summary)) stop("Missing in vitro extra_results seed_summary: ", invitro_summary, call. = FALSE)
  }

  selected <- rbind(
    select_top_seeds(invivo_summary, "in_vivo", top_n = top_n, threshold_N = threshold_N),
    select_top_seeds(invitro_summary, "in_vitro", top_n = top_n, threshold_N = threshold_N)
  )
  utils::write.table(selected, out_selected_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  n_unit <- as_num(cfg[["N_UNIT"]], 22)
  if (!is.finite(n_unit) || n_unit <= 0) n_unit <- 22
  n_dip <- 2 * n_unit
  low_o2_weight <- as_num(argv$joint_composite_low_o2_weight, as_num(cfg[["joint_composite_low_o2_weight"]], 2))
  zero_o2_priority_weight <- as_num(
    argv$joint_composite_zero_o2_priority_weight,
    as_num(cfg[["joint_composite_zero_o2_priority_weight"]], 4)
  )
  auto_settings <- NULL
  if (identical(anchor_mode, "manual")) {
    o2_grid <- split_num(argv$manual_o2_grid, c(0, 0.1, 0.2, 0.5, 1, 2))
    n_grid <- split_num(argv$manual_n_grid, c(44, 66, 88))
    anchor_source <- "manual"
    anchor <- expand.grid(O2 = o2_grid, N = n_grid, KEEP.OUT.ATTRS = FALSE)
    anchor$anchor_source <- anchor_source
    anchor$anchor_role <- "manual_grid"
    anchor$anchor_weight <- 1.0
    anchor$anchor_n_points <- NA_integer_
    anchor$anchor_n_scenarios <- NA_integer_
  } else {
    anchor <- derive_auto_anchor(invivo_run_dir, selected, cfg, argv)
    auto_settings <- attr(anchor, "auto_anchor_settings", exact = TRUE)
    attr(anchor, "auto_anchor_settings") <- NULL
    anchor_source <- "auto_top3_single_fit_trajectory"
    o2_grid <- sort(unique(as.numeric(anchor$O2)))
    n_grid <- sort(unique(as.numeric(anchor$N)))
  }
  anchor <- anchor[order(anchor$O2, anchor$N), , drop = FALSE]
  if (!"anchor_source" %in% names(anchor)) anchor$anchor_source <- anchor_source
  if (!"anchor_role" %in% names(anchor)) anchor$anchor_role <- anchor_source
  if (!"anchor_weight" %in% names(anchor)) anchor$anchor_weight <- 1.0
  anchor$anchor_weight <- suppressWarnings(as.numeric(anchor$anchor_weight))
  anchor$anchor_weight[!is.finite(anchor$anchor_weight) | anchor$anchor_weight <= 0] <- 1.0
  anchor$base_anchor_weight <- anchor_weights(
    anchor,
    low_o2_weight = low_o2_weight,
    zero_o2_priority_weight = zero_o2_priority_weight,
    n_dip = n_dip
  )
  anchor$weight <- anchor$base_anchor_weight * anchor$anchor_weight
  anchor$anchor_mode <- anchor_mode
  anchor$response <- "p_mis_eff(logit);death_load(log10);live_loss_eff(log10)"
  anchor$invivo_run_dir <- invivo_run_dir
  anchor$invitro_run_dir <- invitro_run_dir
  front_cols <- c(
    "O2", "N", "weight", "base_anchor_weight", "anchor_weight", "anchor_mode",
    "anchor_source", "anchor_role", "anchor_n_points", "anchor_n_scenarios",
    "anchor_min_day", "anchor_max_day", "response", "invivo_run_dir", "invitro_run_dir"
  )
  anchor <- anchor[, c(intersect(front_cols, names(anchor)), setdiff(names(anchor), front_cols)), drop = FALSE]
  utils::write.table(anchor, out_anchor_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  cfg[["run_prefix"]] <- as_chr(argv$run_prefix, cfg[["run_prefix"]] %||% basename(out_dir))
  cfg[["append_run_prefix_timestamp"]] <- FALSE
  cfg[["out_root"]] <- normalizePath(as_chr(argv$out_root, dirname(out_dir)), mustWork = FALSE)
  cfg[["data_dir"]] <- normalizePath(as_chr(argv$data_dir, file.path(oxygen_root, "data", "InVivoData_Gemcitabine")), mustWork = FALSE)
  cfg[["joint_composite_penalty"]] <- TRUE
  cfg[["joint_composite_o2_grid"]] <- as.list(as.numeric(o2_grid))
  cfg[["joint_composite_n_grid"]] <- as.numeric(n_grid)
  cfg[["joint_composite_anchor_mode"]] <- "fixed"
  cfg[["joint_composite_anchor_source"]] <- anchor_source
  cfg[["joint_composite_anchor_tsv"]] <- normalizePath(out_anchor_tsv, mustWork = FALSE)
  cfg[["joint_composite_write_bestfit_trajectory_anchor"]] <- FALSE
  cfg[["joint_composite_low_o2_weight"]] <- low_o2_weight
  cfg[["joint_composite_zero_o2_priority_weight"]] <- zero_o2_priority_weight
  cfg[["joint_composite_live_loss_mode"]] <- as_chr(argv$joint_composite_live_loss_mode, as_chr(cfg[["joint_composite_live_loss_mode"]], "transition"))
  if (!is.null(auto_settings)) {
    cfg[["joint_composite_anchor_o2_bin_pct"]] <- auto_settings$o2_bin
    cfg[["joint_composite_anchor_n_bin"]] <- auto_settings$n_bin
    cfg[["joint_composite_anchor_n_quantiles"]] <- as.list(as.numeric(auto_settings$n_quantiles))
    cfg[["joint_composite_anchor_max_dynamic"]] <- auto_settings$max_dynamic
    cfg[["joint_composite_anchor_include_reference"]] <- auto_settings$include_reference
    cfg[["joint_anchor_auto_n_candidates"]] <- auto_settings$n_candidates
    cfg[["joint_anchor_auto_n_anchors"]] <- nrow(anchor)
  }
  cfg[["joint_anchor_invivo_run_dir"]] <- invivo_run_dir
  cfg[["joint_anchor_invitro_run_dir"]] <- invitro_run_dir
  cfg[["joint_anchor_selected_top_n"]] <- top_n

  for (key in c("seeds_file", "parameter_table", "parameters", "init_params_tsv",
                "invitro_parameter_table", "parameter_table_invitro",
                "fit_objects_dir", "flow_density_path")) {
    cfg <- absolutize_if_present(cfg, key, config_dir)
  }
  if (is.null(cfg[["seeds_file"]]) || !nzchar(as_chr(cfg[["seeds_file"]], ""))) {
    cfg[["seeds_file"]] <- file.path(oxygen_root, "data", "O2G_supply_demand", "seeds.csv")
  }

  dir.create(dirname(out_config_yaml), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(cfg, out_config_yaml)

  lines <- c(
    "# Joint Anchor From Single Fits",
    "",
    paste0("- anchor_mode: ", anchor_mode),
    paste0("- anchor_source: ", anchor_source),
    paste0("- invivo_run_dir: ", invivo_run_dir),
    paste0("- invitro_run_dir: ", invitro_run_dir),
    paste0("- top_n: ", top_n),
    paste0("- n_anchor: ", nrow(anchor)),
    paste0("- O2 grid: ", paste(o2_grid, collapse = ", ")),
    paste0("- N grid: ", paste(n_grid, collapse = ", ")),
    if (!is.null(auto_settings)) {
      c(
        paste0("- auto_o2_bin_pct: ", auto_settings$o2_bin),
        paste0("- auto_n_bin: ", auto_settings$n_bin),
        paste0("- auto_n_quantiles: ", paste(auto_settings$n_quantiles, collapse = ", ")),
        paste0("- auto_max_dynamic: ", auto_settings$max_dynamic),
        paste0("- auto_include_reference: ", auto_settings$include_reference),
        paste0("- auto_n_candidates: ", auto_settings$n_candidates)
      )
    } else {
      character()
    },
    paste0("- response: ", anchor$response[[1]]),
    "",
    "## Selected Seeds",
    apply(selected, 1, function(r) {
      paste0("- ", r[["scene"]], " rank ", r[["rank"]], ": ", r[["seed"]], " objective=", signif(as.numeric(r[["objective"]]), 8))
    }),
    "",
    paste0("anchor_tsv: ", out_anchor_tsv),
    paste0("config_yaml: ", out_config_yaml)
  )
  writeLines(lines, out_summary_md)

  cat("Wrote anchor TSV: ", out_anchor_tsv, "\n", sep = "")
  cat("Wrote selected seeds: ", out_selected_tsv, "\n", sep = "")
  cat("Wrote config YAML: ", out_config_yaml, "\n", sep = "")
  cat("Wrote summary: ", out_summary_md, "\n", sep = "")
  invisible(0L)
}

main()
