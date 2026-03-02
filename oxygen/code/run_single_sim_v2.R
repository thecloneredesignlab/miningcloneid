#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Standalone runner for a *single* simulation of the missegregation/WGD model.
#
# Supports two modes:
#   (A) in_vitro  : discrete "passages" (same logic as optimization), measuring
#                  time to grow by POP_GROWTH_FACTOR each passage.
#   (B) in_vivo   : continuous time (days) with optional oxygen schedule and
#                  crowding via run_in_vivo_crowd() (already in model_functions.R).
#
# Typical usage (in vitro):
#   Rscript run_single_sim.R \
#     --mode=in_vitro --O2=0 --init_ploidy=4 --init_layer=post \
#     --R=0.85 --beta=0.35 --pwgd=1e-4 --mr0=0.9 --mr1=0.9 \
#     --pmis1=1e-4 --pmis0=1e-3 --eta=1e-2 \
#     --passages=17 --report_passages=7,17 --pop_growth=10 --dt=0.1 \
#     --out_prefix=sim_4N_O2_0
#
# Typical usage (in vivo):
#   Rscript run_single_sim.R \
#     --mode=in_vivo --init_ploidy=4 --init_layer=post --total_size=1e6 \
#     --R=0.85 --beta=0.35 --pwgd=1e-4 --mr0=0.9 --mr1=0.9 \
#     --pmis1=1e-4 --pmis0=1e-3 --eta=1e-2 \
#     --T_end=28 --sample_days=0,7,14,21,28 \
#     --O2_schedule=0:7:0.2,7:28:1.0 \
#     --K=1e9 --crowding=logistic --out_prefix=invivo_4N_hypoxia7d
#
# Notes
# - Oxygen O2 is scaled 0..1 (matching the core model).
# - For intermediate O2 values, p_mis(O2) is log-linear between pmis0 (O2=0)
#   and pmis1 (O2=1), consistent with plot_misseg_interp().
# -----------------------------------------------------------------------------

# ---- lightweight CLI parsing ----
parse_args <- function(argv) {
  # Supports:
  #   --key=value
  #   --key value
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) {
      stop("Unexpected argument (must start with --): ", a)
    }
    a <- substring(a, 3)
    if (grepl("=", a, fixed = TRUE)) {
      kv <- strsplit(a, "=", fixed = TRUE)[[1]]
      key <- kv[[1]]
      val <- paste(kv[-1], collapse = "=")
      out[[key]] <- val
      i <- i + 1
    } else {
      key <- a
      if (i == length(argv)) stop("Missing value for --", key)
      out[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  out
}

as_num <- function(x) {
  if (is.null(x)) return(NULL)
  as.numeric(x)
}

as_int <- function(x) {
  if (is.null(x)) return(NULL)
  as.integer(x)
}

as_num_vec <- function(x) {
  # Comma-separated numeric list
  if (is.null(x) || nchar(x) == 0) return(NULL)
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

# Parse schedule string like: "0:7:0.2,7:28:1.0" => list(c(t0=0,t1=7,O2=0.2), ...)
parse_O2_schedule <- function(x) {
  if (is.null(x) || nchar(x) == 0) return(list(c(t0=0, t1=Inf, O2=1.0)))
  segs <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- list()
  for (s in segs) {
    parts <- strsplit(s, ":", fixed = TRUE)[[1]]
    if (length(parts) != 3) stop("Bad O2_schedule segment (need t0:t1:O2): ", s)
    out[[length(out) + 1]] <- c(t0 = as.numeric(parts[[1]]),
                               t1 = as.numeric(parts[[2]]),
                               O2 = as.numeric(parts[[3]]))
  }
  out
}

# Find script dir so relative paths (like code/model_functions.R) work.
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grep("^--file=", cmd)]
  if (length(file_arg) == 0) return(getwd())
  this_file <- sub("^--file=", "", file_arg[[1]])
  dirname(normalizePath(this_file))
}

# ---- single simulation (in vitro passages) ----
run_one_sim_in_vitro <- function(run_params,
                                O2,
                                init_state,
                                N_UNIT = 22L,
                                N_MIN = 22L,
                                N_MAX = 154L,
                                DT = 0.1,
                                POP_GROWTH_FACTOR = 10.0,
                                PASSAGES_TO_RUN = 17L,
                                REPORT_PASSAGES = c(7L, 17L),
                                boundary = "absorb_minmax") {

  grid_pre  <- N_MIN:N_MAX
  grid_post <- N_MIN:N_MAX
  R0 <- length(grid_pre)
  R1 <- length(grid_post)

  # Allow O2 in [0,1] and interpolate p_mis log-linearly between fitted endpoints
  O2 <- max(0, min(1, O2))
  logp <- (1 - O2) * log10(run_params$pmis_O2_0) + O2 * log10(run_params$pmis_O2_1)
  pmis_const <- 10^logp

  # Build generator G
  lambda0_vec <- growth_lambda(O2, grid_pre,  R = run_params$R, beta = run_params$beta,
                               eta = run_params$eta, N_unit = N_UNIT)
  lambda1_vec <- growth_lambda(O2, grid_post, R = run_params$R, beta = run_params$beta,
                               eta = run_params$eta, N_unit = N_UNIT)

  G <- .build_G_with_WGD(
    N0min = N_MIN, N0max = N_MAX,
    lambda0_vec = lambda0_vec,
    p0_vec = pmis_const,
    wgd_prob_vec = run_params$pwgd,
    N1min = N_MIN, N1max = N_MAX,
    lambda1_vec = lambda1_vec,
    p1_vec = pmis_const,
    mr_lethality0 = run_params$mr_lethality0,
    mr_lethality1 = run_params$mr_lethality1,
    boundary = boundary
  )

  x_current <- as.numeric(init_state)
  if (length(x_current) != (R0 + R1)) {
    stop("init_state has length ", length(x_current), " but expected ", (R0 + R1))
  }
  if (sum(x_current) <= 0) stop("init_state has zero mass")

  all_results_list <- list()
  sim_passage_times <- numeric(PASSAGES_TO_RUN)

  for (p in seq_len(PASSAGES_TO_RUN)) {
    pop_start  <- sum(x_current)
    pop_target <- pop_start * POP_GROWTH_FACTOR
    time_in_passage <- 0.0

    while (sum(x_current) < pop_target) {
      x_current <- step_dt(G, x_current, DT, 1L)
      time_in_passage <- time_in_passage + DT

      # safety breaks
      if (sum(x_current) < pop_start * 1e-12) break
      if (time_in_passage > 1e4) break
    }

    sim_passage_times[p] <- time_in_passage

    if (p %in% REPORT_PASSAGES) {
      pop_total <- sum(x_current)
      dist_df <- data.frame(
        passage = p,
        layer = c(rep("pre", R0), rep("post", R1)),
        N = c(grid_pre, grid_post),
        fraction = x_current / pop_total,
        pop = pop_total
      )
      all_results_list[[length(all_results_list) + 1]] <- dist_df
    }

    # dilution/back-to-start size (same as optimization)
    x_current <- x_current / sum(x_current) * pop_start
  }

  all_dists <- do.call(rbind, all_results_list)
  all_passage_times <- data.frame(passage = seq_len(PASSAGES_TO_RUN), duration = sim_passage_times)

  list(all_dists = all_dists, all_passage_times = all_passage_times)
}

# ---- summary helpers ----
collapse_layers_to_ploidy <- function(dists_df, N_UNIT = 22L) {
  # Aggregate pre+post by N, convert to P = N/N_UNIT.
  # Works for both in_vitro ("passage" column) and in_vivo ("day" column).
  if ("passage" %in% names(dists_df)) {
    tcol <- "passage"
  } else if ("day" %in% names(dists_df)) {
    tcol <- "day"
  } else {
    stop("dists_df must contain either a 'passage' column (in vitro) or a 'day' column (in vivo).")
  }

  tmp <- dplyr::as_tibble(dists_df)
  tmp$time <- tmp[[tcol]]
  tmp %>%
  dplyr::group_by(time, N) %>%
  dplyr::summarise(fraction = sum(fraction), .groups = "drop") %>%
  dplyr::mutate(P = N / N_UNIT)
}


# ---- weighted quantiles and ploidy timecourse summaries ----
weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  stopifnot(length(x) == length(w))
  ord <- order(x)
  x <- x[ord]; w <- w[ord]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) {
    x[which(cw >= p)[1]]
    })
}

summarize_ploidy_timecourse <- function(dP) {
  # dP must have columns: time, P, fraction
  dplyr::as_tibble(dP) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      mean_ploidy   = sum(P * fraction),
      q25_ploidy    = weighted_quantile(P, fraction, probs = 0.25)[1],
      median_ploidy = weighted_quantile(P, fraction, probs = 0.50)[1],
      q75_ploidy    = weighted_quantile(P, fraction, probs = 0.75)[1],
      frac_gt_2_5   = sum(fraction[P > 2.5]),
      frac_gt_3_0   = sum(fraction[P > 3.0]),
      .groups = "drop"
    )
}


# ---- main ----
argv <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(argv)

mode <- if (!is.null(opt$mode)) opt$mode else "in_vitro"

# locate and source model_functions.R
script_dir <- get_script_dir()
model_path <- if (!is.null(opt$model_functions)) opt$model_functions else file.path(script_dir, "scr", "model_functions_ploidy_buffer.R")
if (!file.exists(model_path)) {
  # try common repo layout
  alt <- file.path(script_dir, "code", "scr", "model_functions_ploidy_buffer.R")
  if (file.exists(alt)) model_path <- alt
}
if (!file.exists(model_path)) stop("Could not find model_functions.R at: ", model_path)
source(model_path)

# ---- parameters (real scale) ----
# Required params: R, beta, pwgd, mr0, mr1, pmis1, pmis0, eta
# Aliases: mr0/mr1 => mr_lethality0/1; pmis1/pmis0 => pmis_O2_1/0
required <- c("R", "beta", "pwgd", "mr0", "mr1", "pmis1", "pmis0", "eta")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop(
    "Missing required args: ", paste(missing, collapse = ", "), "\n",
    "Example: --R=0.85 --beta=0.35 --pwgd=1e-4 --mr0=0.9 --mr1=0.9 --pmis1=1e-4 --pmis0=1e-3 --eta=1e-2"
  )
}

run_params <- list(
  R             = as.numeric(opt$R),
  beta          = as.numeric(opt$beta),
  pwgd          = as.numeric(opt$pwgd),
  mr_lethality0 = as.numeric(opt$mr0),
  mr_lethality1 = as.numeric(opt$mr1),
  pmis_O2_1     = as.numeric(opt$pmis1),
  pmis_O2_0     = as.numeric(opt$pmis0),
  eta           = as.numeric(opt$eta)
)

# ---- common simulation options ----
N_UNIT <- if (!is.null(opt$N_UNIT)) as.integer(opt$N_UNIT) else 22L
N_MIN  <- if (!is.null(opt$N_MIN))  as.integer(opt$N_MIN)  else 22L
N_MAX  <- if (!is.null(opt$N_MAX))  as.integer(opt$N_MAX)  else 154L
DT     <- if (!is.null(opt$dt))     as.numeric(opt$dt)     else 0.1

grid_pre  <- N_MIN:N_MAX
grid_post <- N_MIN:N_MAX

# ---- init state ----
init_ploidy <- if (!is.null(opt$init_ploidy)) as.integer(opt$init_ploidy) else 2L
init_layer  <- if (!is.null(opt$init_layer)) opt$init_layer else "pre"
init_layer  <- match.arg(init_layer, c("pre", "post"))

total_size <- if (!is.null(opt$total_size)) as.numeric(opt$total_size) else 1.0

init_state <- make_init_state(
  grid_pre = grid_pre,
  grid_post = grid_post,
  ploidy = init_ploidy,
  layer  = init_layer,
  N_UNIT = N_UNIT,
  total_size = total_size
)

# Optionally replace delta-init with empirical init from a ploidy_distribution.csv
# expected columns: id, ploidy (numeric), ...
# filter by: passage==0 and ploidy label (2N/4N) like in run_optim.R
if (!is.null(opt$init_ploidy_csv)) {
  suppressPackageStartupMessages(library(data.table))
  x <- data.table::fread(opt$init_ploidy_csv)
  if (!all(c("id", "ploidy") %in% colnames(x))) {
    stop("init_ploidy_csv must include columns: id, ploidy")
  }
  x$P <- x$ploidy
  x$hypoxia <- grepl("_O", x$id)
  x$ploidy_label <- "2N"
  x$ploidy_label[grepl("4N", x$id)] <- "4N"
  x$passage <- 0
  x$passage[grepl("_A7", x$id) & x$hypoxia] <- 7
  x$passage[grepl("_A17", x$id) | grepl("_A19", x$id)] <- 17

  init_label <- if (init_ploidy == 4) "4N" else "2N"
  init_P_values <- x$P[x$passage == 0 & x$ploidy_label == init_label]
  if (length(init_P_values) == 0) stop("No empirical init cells found for ", init_label, " at passage 0")

  # build fractions
  x_pre  <- rep(0, length(grid_pre));  names(x_pre)  <- grid_pre
  x_post <- rep(0, length(grid_post)); names(x_post) <- grid_post
  if (init_layer == "pre") {
    x_pre <- create_initial_dist(init_P_values, grid_pre, N_UNIT)
  } else {
    x_post <- create_initial_dist(init_P_values, grid_post, N_UNIT)
  }
  init_state <- c(x_pre, x_post)
  init_state <- init_state / sum(init_state) * total_size
}

# ---- output prefix ----
out_prefix <- if (!is.null(opt$out_prefix)) opt$out_prefix else "single_sim"

# ---- run ----
if (mode == "in_vitro") {

  O2 <- if (!is.null(opt$O2)) as.numeric(opt$O2) else 1.0
  PASSAGES <- if (!is.null(opt$passages)) as.integer(opt$passages) else 17L
  POP_GROWTH <- if (!is.null(opt$pop_growth)) as.numeric(opt$pop_growth) else 10.0
  REPORT_PASSAGES <- as.integer(as_num_vec(opt$report_passages))
  if (is.null(REPORT_PASSAGES)) REPORT_PASSAGES <- c(7L, PASSAGES)

  sim <- run_one_sim_in_vitro(
    run_params = run_params,
    O2 = O2,
    init_state = init_state,
    N_UNIT = N_UNIT,
    N_MIN = N_MIN,
    N_MAX = N_MAX,
    DT = DT,
    POP_GROWTH_FACTOR = POP_GROWTH,
    PASSAGES_TO_RUN = PASSAGES,
    REPORT_PASSAGES = REPORT_PASSAGES
  )

  # write outputs
  utils::write.csv(sim$all_dists, paste0(out_prefix, "_dists.csv"), row.names = FALSE)
  utils::write.csv(sim$all_passage_times, paste0(out_prefix, "_passage_times.csv"), row.names = FALSE)

  # quick plot if ggplot2 is available (it is loaded by model_functions.R)
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  dP <- collapse_layers_to_ploidy(sim$all_dists, N_UNIT)
  # plot final reported passage by default
  final_p <- max(dP$time)
  p_hist <- ggplot(dP %>% filter(time == final_p), aes(x = P, weight = fraction)) +
    geom_histogram(bins = 60, color = "black", fill = "grey80") +
    theme_bw() +
    labs(title = paste0("Ploidy distribution (passage ", final_p, ")"), x = "Ploidy (N / N_UNIT)", y = "Probability mass")
  ggsave(filename = paste0(out_prefix, "_ploidy_hist.png"), plot = p_hist, width = 7, height = 4, dpi = 200)

  # ploidy over time (passage-indexed, with optional mapping to cumulative time)
  pt <- summarize_ploidy_timecourse(dP)
  utils::write.csv(pt, paste0(out_prefix, "_ploidy_timecourse.csv"), row.names = FALSE)

  # map passage -> cumulative time (sum of passage durations) for an actual time axis
  pt_time <- sim$all_passage_times %>%
    dplyr::group_by(passage) %>%
    dplyr::summarise(duration = mean(duration), .groups="drop") %>%
    dplyr::arrange(passage) %>%
    dplyr::mutate(cum_time = cumsum(duration))

  pt_plotdf <- pt %>%
    dplyr::left_join(pt_time, by = c("time" = "passage")) %>%
    dplyr::mutate(x_time = dplyr::if_else(is.na(cum_time), as.numeric(time), cum_time))

  p_ploidy <- ggplot(pt_plotdf, aes(x = x_time, y = mean_ploidy)) +
    geom_ribbon(aes(ymin = q25_ploidy, ymax = q75_ploidy), fill = "grey85") +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 2) +
    theme_bw() +
    labs(
      title = "Ploidy over time",
      x = "Time (cumulative passage time; falls back to passage index if missing)",
      y = "Ploidy (N / N_UNIT)"
    )
  ggsave(filename = paste0(out_prefix, "_ploidy_timecourse.png"), plot = p_ploidy, width = 7, height = 4, dpi = 200)

  cat("\nWrote:\n",
      " - ", paste0(out_prefix, "_dists.csv"), "\n",
      " - ", paste0(out_prefix, "_passage_times.csv"), "\n",
      " - ", paste0(out_prefix, "_ploidy_hist.png"), "\n")

} else if (mode == "in_vivo") {

  T_end <- if (!is.null(opt$T_end)) as.numeric(opt$T_end) else 28
  sample_days <- as_num_vec(opt$sample_days)
  if (is.null(sample_days)) sample_days <- c(0, 7, 14, 21, 28)

  O2_schedule <- parse_O2_schedule(opt$O2_schedule)

  K <- if (!is.null(opt$K)) as.numeric(opt$K) else 1e9
  crowding <- if (!is.null(opt$crowding)) opt$crowding else "logistic"

  sim <- run_in_vivo_crowd(
    run_params = run_params,
    O2_schedule = O2_schedule,
    T_end = T_end,
    sample_days = sample_days,
    N_UNIT = N_UNIT,
    DT = DT,
    K = K,
    crowding = crowding,
    grid_pre = grid_pre,
    grid_post = grid_post,
    init_state = init_state
  )

  utils::write.csv(sim$all_dists,  paste0(out_prefix, "_dists.csv"), row.names = FALSE)
  utils::write.csv(sim$tumor_size, paste0(out_prefix, "_tumor_size.csv"), row.names = FALSE)

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  # tumor size trace plot
  p_size <- ggplot(sim$tumor_size, aes(x = day, y = Ntot)) +
    geom_line(color = "black", linewidth = 1) +
    theme_bw() +
    labs(title = "Tumor size (total population)", x = "Day", y = "Total population")
  ggsave(filename = paste0(out_prefix, "_tumor_size.png"), plot = p_size, width = 7, height = 4, dpi = 200)

  # ploidy distribution at last sample day
  dP <- collapse_layers_to_ploidy(sim$all_dists, N_UNIT)
  final_day <- max(dP$time)
  p_hist <- ggplot(dP %>% filter(time == final_day), aes(x = P, weight = fraction)) +
    geom_histogram(bins = 60, color = "black", fill = "grey80") +
    theme_bw() +
    labs(title = paste0("Ploidy distribution (day ", final_day, ")"), x = "Ploidy (N / N_UNIT)", y = "Probability mass")
  ggsave(filename = paste0(out_prefix, "_ploidy_hist.png"), plot = p_hist, width = 7, height = 4, dpi = 200)

  # ploidy over time (sample days)
  pt <- summarize_ploidy_timecourse(dP)
  utils::write.csv(pt, paste0(out_prefix, "_ploidy_timecourse.csv"), row.names = FALSE)

  p_ploidy <- ggplot(pt, aes(x = time, y = mean_ploidy)) +
    geom_ribbon(aes(ymin = q25_ploidy, ymax = q75_ploidy), fill = "grey85") +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 2) +
    theme_bw() +
    labs(title = "Ploidy over time", x = "Day", y = "Ploidy (N / N_UNIT)")
  ggsave(filename = paste0(out_prefix, "_ploidy_timecourse.png"), plot = p_ploidy, width = 7, height = 4, dpi = 200)

  # optional combined (tumor size + mean ploidy) in one PNG for quick inspection
  png(filename = paste0(out_prefix, "_trajectory_and_ploidy.png"), width = 900, height = 700)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  plot(sim$tumor_size$day, sim$tumor_size$Ntot, type = "l", lwd = 2,
       xlab = "Day", ylab = "Total population", main = "Tumor burden (total population)")
  plot(pt$time, pt$mean_ploidy, type = "b", pch = 16, lwd = 2,
       xlab = "Day", ylab = "Mean ploidy (N/N_UNIT)", main = "Mean ploidy over time")
  dev.off()

  cat("\nWrote:\n",
      " - ", paste0(out_prefix, "_dists.csv"), "\n",
      " - ", paste0(out_prefix, "_tumor_size.csv"), "\n",
      " - ", paste0(out_prefix, "_tumor_size.png"), "\n",
      " - ", paste0(out_prefix, "_ploidy_hist.png"), "\n")

} else {
  stop("Unknown mode: ", mode, " (use in_vitro or in_vivo)")
}
