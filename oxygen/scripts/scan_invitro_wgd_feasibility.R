#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  out
}

parse_num_grid <- function(value, default) {
  if (is.null(value) || !nzchar(value)) return(default)
  ans <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (!length(ans) || any(!is.finite(ans))) stop("Grid values must be finite numbers.")
  unique(ans)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- normalizePath(if (is.null(args$repo_root)) getwd() else args$repo_root, mustWork = TRUE)
oxygen_root <- normalizePath(file.path(repo_root, "oxygen"), mustWork = TRUE)
workflow_root <- normalizePath(file.path(oxygen_root, "code", "O2_supply_demand_MAP"), mustWork = TRUE)

source(file.path(workflow_root, "util", "o2_supply_demand_map_shared.R"), local = environment())
source(file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R"), local = environment())
source(file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R"), local = environment())
source(
  file.path(workflow_root, "util", "o2_supply_demand_map_invitro_utils.R"),
  local = environment(),
  chdir = TRUE
)

fit_result_path <- normalizePath(
  if (is.null(args$fit_result)) stop("--fit_result is required") else args$fit_result,
  mustWork = TRUE
)
out_dir <- normalizePath(
  if (is.null(args$out_dir)) stop("--out_dir is required") else args$out_dir,
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_wgd_values <- parse_num_grid(
  args$p_wgd_values,
  c(1e-8, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
)
p_misseg_values <- parse_num_grid(
  args$p_misseg_values,
  c(0.0025, 0.005, 0.01, 0.02, 0.04, 0.06, 0.0803768416171454)
)
k_o_mis_values <- parse_num_grid(args$k_o_mis_values, base::numeric())
buffer_smax_values <- parse_num_grid(args$buffer_smax_values, base::numeric())
buffer_beta_values <- parse_num_grid(args$buffer_beta_values, base::numeric())
buffer_n_exp_values <- parse_num_grid(args$buffer_n_exp_values, base::numeric())
n_cores <- max(1L, as.integer(if (is.null(args$n_cores)) 1L else args$n_cores))

fit_result <- readRDS(fit_result_path)
cfg <- fit_result$cfg
base_params <- fit_result$best_params
if (is.null(cfg) || is.null(base_params)) stop("fit_result must contain cfg and best_params")
if (!length(k_o_mis_values)) k_o_mis_values <- as.numeric(base_params$k_o_mis)
if (!length(buffer_smax_values)) buffer_smax_values <- as.numeric(base_params$buffer_smax)
if (!length(buffer_beta_values)) buffer_beta_values <- as.numeric(base_params$buffer_beta)
if (!length(buffer_n_exp_values)) buffer_n_exp_values <- as.numeric(base_params$buffer_n_exp)

fit_objects <- ivt_load_fit_objects(
  repo_root = oxygen_root,
  fit_objects_dir = normalizePath(file.path(oxygen_root, "ploidyOxygen", "data", "fit_objects"), mustWork = TRUE),
  flow_csv_path = normalizePath(file.path(oxygen_root, "data", "g0g1_ploidy_density_grid.csv"), mustWork = TRUE),
  death_data_path = normalizePath(
    file.path(repo_root, "data", "InVitroData", "sum159_dead_cell_endpoint_likelihood_ready_20260731.tsv"),
    mustWork = TRUE
  )
)

regular_mass <- function(x, y, lower = -Inf, upper = Inf) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- pmax(y[keep], 0)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  dx <- stats::median(diff(x))
  total <- sum(y) * dx
  if (!is.finite(total) || total <= 0) return(NA_real_)
  sum(y[x >= lower & x < upper]) * dx / total
}

flow_metrics_one <- function(df) {
  df <- df[order(df$ploidy), , drop = FALSE]
  x <- as.numeric(df$ploidy)
  y <- pmax(as.numeric(df$density), 0)
  dx <- stats::median(diff(x))
  total <- sum(y) * dx
  if (!is.finite(total) || total <= 0) {
    return(data.frame(
      mean_ploidy = NA_real_, sd_ploidy = NA_real_, bridge_mass = NA_real_,
      wgd_mass = NA_real_, mass_3_5 = NA_real_, tail_gt4_5 = NA_real_,
      tail_gt5 = NA_real_, low_peak = NA_real_,
      high_peak = NA_real_, peak_ratio = NA_real_, valley_ratio = NA_real_,
      has_target_bimodality = FALSE
    ))
  }
  y <- y / total
  mean_x <- sum(x * y) * dx
  sd_x <- sqrt(max(0, sum((x - mean_x)^2 * y) * dx))
  local_idx <- which(
    y[2:(length(y) - 1L)] > y[1:(length(y) - 2L)] &
      y[2:(length(y) - 1L)] >= y[3:length(y)]
  ) + 1L
  low_idx <- local_idx[x[local_idx] >= 1.4 & x[local_idx] <= 2.8]
  high_idx <- local_idx[x[local_idx] >= 3.3 & x[local_idx] <= 5.3]
  low_idx <- if (length(low_idx)) low_idx[[which.max(y[low_idx])]] else NA_integer_
  high_idx <- if (length(high_idx)) high_idx[[which.max(y[high_idx])]] else NA_integer_
  low_peak <- if (is.finite(low_idx)) x[[low_idx]] else NA_real_
  high_peak <- if (is.finite(high_idx)) x[[high_idx]] else NA_real_
  peak_ratio <- if (is.finite(low_peak) && low_peak > 0 && is.finite(high_peak)) high_peak / low_peak else NA_real_
  valley_ratio <- NA_real_
  high_relative_height <- NA_real_
  if (is.finite(low_idx) && is.finite(high_idx) && high_idx > low_idx) {
    valley_idx <- seq.int(low_idx, high_idx)
    valley_ratio <- min(y[valley_idx]) / min(y[[low_idx]], y[[high_idx]])
    high_relative_height <- y[[high_idx]] / max(y)
  }
  target_bimodal <- is.finite(peak_ratio) && peak_ratio >= 1.75 && peak_ratio <= 2.25 &&
    is.finite(valley_ratio) && valley_ratio <= 0.5 &&
    is.finite(high_relative_height) && high_relative_height >= 0.08
  data.frame(
    mean_ploidy = mean_x,
    sd_ploidy = sd_x,
    bridge_mass = regular_mass(x, y, 2.5, 3.5),
    wgd_mass = regular_mass(x, y, 3.5, Inf),
    mass_3_5 = regular_mass(x, y, 3, 5),
    tail_gt4_5 = regular_mass(x, y, 4.5, Inf),
    tail_gt5 = regular_mass(x, y, 5, Inf),
    low_peak = low_peak,
    high_peak = high_peak,
    peak_ratio = peak_ratio,
    valley_ratio = valley_ratio,
    has_target_bimodality = target_bimodal
  )
}

summarize_flow <- function(overlay, grid_id, p_wgd, p_misseg) {
  pred <- overlay[overlay$series == "Predicted", , drop = FALSE]
  ids <- unique(as.character(pred$segment_id))
  bind_rows(lapply(ids, function(id) {
    z <- pred[pred$segment_id == id, , drop = FALSE]
    cbind(
      data.frame(
        grid_id = grid_id,
        p_wgd = p_wgd,
        p_misseg = p_misseg,
        segment_id = id,
        cohort = as.character(z$cohort[[1L]]),
        lineage_id = as.character(z$lineage_id[[1L]]),
        oxygen_pct = as.numeric(z$oxygen_pct[[1L]]),
        stringsAsFactors = FALSE
      ),
      flow_metrics_one(z)
    )
  }))
}

grid <- expand.grid(
  p_wgd = p_wgd_values,
  p_misseg = p_misseg_values,
  k_o_mis = k_o_mis_values,
  buffer_smax = buffer_smax_values,
  buffer_beta = buffer_beta_values,
  buffer_n_exp = buffer_n_exp_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$grid_id <- sprintf("grid_%03d", seq_len(nrow(grid)))

evaluate_one <- function(i) {
  row <- grid[i, , drop = FALSE]
  run_params <- base_params
  run_params$p_wgd <- row$p_wgd[[1L]]
  run_params$p_misseg <- row$p_misseg[[1L]]
  run_params$k_o_mis <- row$k_o_mis[[1L]]
  run_params$buffer_smax <- row$buffer_smax[[1L]]
  run_params$buffer_beta <- row$buffer_beta[[1L]]
  run_params$buffer_n_exp <- row$buffer_n_exp[[1L]]
  comp <- tryCatch(
    ivt_objective_components(
      run_params = run_params,
      fit_objects = fit_objects,
      cfg = cfg,
      fallback_max_passage_days = 14,
      growth_weight = 1,
      ploidy_weight = as.numeric(fit_result$ploidy_weight %||% 1),
      flow_weight = 1,
      death_weight = as.numeric(fit_result$death_weight %||% 1),
      passage_time_weight = as.numeric(fit_result$passage_time_weight %||% 0.25),
      passage_time_tolerance_days = as.numeric(fit_result$passage_time_tolerance_days %||% 1),
      passage_time_sigma_days = as.numeric(fit_result$passage_time_sigma_days %||% 1),
      passage_time_df = as.numeric(fit_result$passage_time_df %||% 4),
      sigma_death_logit = as.numeric(fit_result$sigma_death_logit %||% 0.75),
      death_fraction_eps = as.numeric(fit_result$death_fraction_eps %||% 1e-4)
    ),
    error = function(e) e
  )
  if (inherits(comp, "error")) {
    return(list(
      summary = data.frame(
        grid_id = row$grid_id, p_wgd = row$p_wgd, p_misseg = row$p_misseg,
        k_o_mis = row$k_o_mis,
        buffer_smax = row$buffer_smax, buffer_beta = row$buffer_beta,
        buffer_n_exp = row$buffer_n_exp,
        status = "ERROR", error = conditionMessage(comp), objective = NA_real_,
        growth_loglik = NA_real_, ploidy_loglik = NA_real_, flow_loglik = NA_real_,
        death_loglik = NA_real_, passage_time_loglik = NA_real_,
        protocol_feasibility_status = NA_character_, n_insufficient_boundaries = NA_integer_
      ),
      flow = data.frame()
    ))
  }
  list(
    summary = data.frame(
      grid_id = row$grid_id, p_wgd = row$p_wgd, p_misseg = row$p_misseg,
      k_o_mis = row$k_o_mis,
      buffer_smax = row$buffer_smax, buffer_beta = row$buffer_beta,
      buffer_n_exp = row$buffer_n_exp,
      status = "OK", error = "", objective = comp$objective,
      growth_loglik = comp$growth_loglik, ploidy_loglik = comp$ploidy_loglik,
      flow_loglik = comp$flow_loglik, death_loglik = comp$death_loglik,
      passage_time_loglik = comp$passage_time_loglik,
      protocol_feasibility_status = comp$protocol_feasibility_status,
      n_insufficient_boundaries = comp$n_insufficient_boundaries
    ),
    flow = summarize_flow(comp$flow_overlay_df, row$grid_id, row$p_wgd, row$p_misseg)
  )
}

indices <- seq_len(nrow(grid))
rows <- if (n_cores > 1L && .Platform$OS.type == "unix") {
  parallel::mclapply(indices, evaluate_one, mc.cores = n_cores, mc.preschedule = FALSE)
} else {
  lapply(indices, evaluate_one)
}

summary_df <- bind_rows(lapply(rows, `[[`, "summary"))
flow_df <- bind_rows(lapply(rows, `[[`, "flow"))

aggregate_one <- function(grid_id) {
  z <- flow_df[flow_df$grid_id == grid_id, , drop = FALSE]
  hyp2 <- z[grepl("^2N-O[12]-A(12|18|23)$", z$segment_id), , drop = FALSE]
  hyp2_terminal <- z[grepl("^2N-O[12]-A23$", z$segment_id), , drop = FALSE]
  c2 <- z[z$segment_id == "2N-C-A12", , drop = FALSE]
  h4a22 <- z[grepl("^4N-O[12]-A22$", z$segment_id), , drop = FALSE]
  h4a1 <- z[grepl("^4N-O[12]-A6$", z$segment_id), , drop = FALSE]
  data.frame(
    grid_id = grid_id,
    n_2N_hypoxia_target_bimodal = sum(hyp2$has_target_bimodality, na.rm = TRUE),
    n_2N_terminal_target_bimodal = sum(hyp2_terminal$has_target_bimodality, na.rm = TRUE),
    mean_2N_hypoxia_bridge_mass = mean(hyp2$bridge_mass, na.rm = TRUE),
    mean_2N_hypoxia_wgd_mass = mean(hyp2$wgd_mass, na.rm = TRUE),
    mean_2N_hypoxia_valley_ratio = mean(hyp2$valley_ratio, na.rm = TRUE),
    mean_2N_terminal_bridge_mass = mean(hyp2_terminal$bridge_mass, na.rm = TRUE),
    mean_2N_terminal_wgd_mass = mean(hyp2_terminal$wgd_mass, na.rm = TRUE),
    mean_2N_terminal_tail_gt4_5 = mean(hyp2_terminal$tail_gt4_5, na.rm = TRUE),
    mean_2N_control_wgd_mass = mean(c2$wgd_mass, na.rm = TRUE),
    mean_4N_A22_ploidy = mean(h4a22$mean_ploidy, na.rm = TRUE),
    mean_4N_A22_sd = mean(h4a22$sd_ploidy, na.rm = TRUE),
    mean_4N_A22_mass_3_5 = mean(h4a22$mass_3_5, na.rm = TRUE),
    mean_4N_A22_tail_gt5 = mean(h4a22$tail_gt5, na.rm = TRUE),
    delta_4N_mean_A22_vs_A6 = mean(h4a22$mean_ploidy, na.rm = TRUE) - mean(h4a1$mean_ploidy, na.rm = TRUE)
  )
}

agg_df <- if (nrow(flow_df)) {
  bind_rows(lapply(unique(flow_df$grid_id), aggregate_one))
} else {
  data.frame(
    grid_id = character(),
    n_2N_hypoxia_target_bimodal = integer(),
    n_2N_terminal_target_bimodal = integer(),
    mean_2N_hypoxia_bridge_mass = numeric(),
    mean_2N_hypoxia_wgd_mass = numeric(),
    mean_2N_hypoxia_valley_ratio = numeric(),
    mean_2N_terminal_bridge_mass = numeric(),
    mean_2N_terminal_wgd_mass = numeric(),
    mean_2N_terminal_tail_gt4_5 = numeric(),
    mean_2N_control_wgd_mass = numeric(),
    mean_4N_A22_ploidy = numeric(),
    mean_4N_A22_sd = numeric(),
    mean_4N_A22_mass_3_5 = numeric(),
    mean_4N_A22_tail_gt5 = numeric(),
    delta_4N_mean_A22_vs_A6 = numeric()
  )
}
summary_df <- left_join(summary_df, agg_df, by = "grid_id")
summary_df$passes_shape_screen <- with(
  summary_df,
  status == "OK" & protocol_feasibility_status == "PASS" &
    n_2N_hypoxia_target_bimodal >= 4 &
    n_2N_terminal_target_bimodal == 2 &
    mean_2N_terminal_bridge_mass <= 0.10 &
    mean_2N_terminal_wgd_mass >= 0.15 &
    mean_2N_terminal_tail_gt4_5 <= 0.10 &
    mean_2N_control_wgd_mass <= 0.02 &
    mean_4N_A22_mass_3_5 >= 0.90 &
    mean_4N_A22_tail_gt5 <= 0.02 &
    delta_4N_mean_A22_vs_A6 <= 0.25
)
summary_df <- summary_df[order(
  !summary_df$passes_shape_screen,
  -summary_df$n_2N_hypoxia_target_bimodal,
  summary_df$mean_2N_hypoxia_bridge_mass,
  summary_df$objective
), , drop = FALSE]

write.table(summary_df, file.path(out_dir, "wgd_misseg_feasibility_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(flow_df, file.path(out_dir, "wgd_misseg_flow_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(
  data.frame(
    key = c(
      "fit_result", "p_wgd_values", "p_misseg_values", "k_o_mis_values", "buffer_smax_values",
      "buffer_beta_values", "buffer_n_exp_values", "n_cores", "grid_size", "generated_at"
    ),
    value = c(
      fit_result_path, paste(p_wgd_values, collapse = ","), paste(p_misseg_values, collapse = ","),
      paste(k_o_mis_values, collapse = ","),
      paste(buffer_smax_values, collapse = ","), paste(buffer_beta_values, collapse = ","),
      paste(buffer_n_exp_values, collapse = ","), n_cores, nrow(grid),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  ),
  file.path(out_dir, "wgd_misseg_feasibility_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("WGD feasibility scan written to: ", normalizePath(out_dir, mustWork = TRUE))
