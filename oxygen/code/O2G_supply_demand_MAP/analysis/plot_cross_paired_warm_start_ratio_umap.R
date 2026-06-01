#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(uwot)
})

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
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

as_chr <- function(x, default = "") {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) default else as.character(x[[1]])
}

as_int <- function(x, default) {
  val <- suppressWarnings(as.integer(x[[1]]))
  if (length(val) && is.finite(val)) val else default
}

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]), mustWork = FALSE))
project_root <- normalizePath(file.path(script_dir, "../../../.."), mustWork = FALSE)

PARAMETERS <- c(
  "O2_crit",
  "alpha_o2",
  "mu_hp",
  "p_misseg",
  "k_o_mis",
  "buffer_smax",
  "buffer_beta",
  "buffer_n_exp",
  "n_O",
  "gamma_growth",
  "lam_max",
  "p_mis_base",
  "p_wgd",
  "gamma_mu"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))
results_dir <- normalizePath(
  as_chr(args$results_dir, file.path(project_root, "oxygen/results/warm_start_from_separate_500seed")),
  mustWork = FALSE
)
best_result_dir <- normalizePath(
  as_chr(args$best_result_dir, file.path(results_dir, "best_result")),
  mustWork = FALSE
)
out_dir <- normalizePath(
  as_chr(args$out_dir, file.path(results_dir, "cross_paired_ratio_umap")),
  mustWork = FALSE
)
top_n <- as_int(args$top_n, 10L)
umap_seed <- as_int(args$umap_seed, 1L)

invivo_csv <- file.path(best_result_dir, "invivo_objective_ranking.csv")
invitro_csv <- file.path(best_result_dir, "invitro_objective_ranking.csv")
if (!file.exists(invivo_csv)) stop("Missing invivo ranking CSV: ", invivo_csv, call. = FALSE)
if (!file.exists(invitro_csv)) stop("Missing invitro ranking CSV: ", invitro_csv, call. = FALSE)

invivo <- read_csv(invivo_csv, show_col_types = FALSE)
invitro <- read_csv(invitro_csv, show_col_types = FALSE)

missing_vivo <- setdiff(PARAMETERS, names(invivo))
missing_vitro <- setdiff(PARAMETERS, names(invitro))
if (length(missing_vivo)) stop("Missing in vivo parameter columns: ", paste(missing_vivo, collapse = ", "), call. = FALSE)
if (length(missing_vitro)) stop("Missing in vitro parameter columns: ", paste(missing_vitro, collapse = ", "), call. = FALSE)

invivo_top <- invivo[order(invivo$rank), ][seq_len(min(top_n, nrow(invivo))), , drop = FALSE]
invitro_top <- invitro[order(invitro$rank), ][seq_len(min(top_n, nrow(invitro))), , drop = FALSE]

pair_rows <- list()
idx <- 1L
for (i in seq_len(nrow(invivo_top))) {
  for (j in seq_len(nrow(invitro_top))) {
    vivo_vals <- as.numeric(invivo_top[i, PARAMETERS])
    vitro_vals <- as.numeric(invitro_top[j, PARAMETERS])
    ratio_vals <- vivo_vals / vitro_vals
    log_ratio_vals <- log10(ratio_vals)
    if (any(!is.finite(log_ratio_vals))) {
      bad <- PARAMETERS[!is.finite(log_ratio_vals)]
      stop(
        "Non-finite log10 ratio for pair V", invivo_top$rank[[i]], "-I", invitro_top$rank[[j]],
        ": ", paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
    row <- data.frame(
      pair_id = sprintf("V%02d-I%02d", invivo_top$rank[[i]], invitro_top$rank[[j]]),
      invivo_rank = as.integer(invivo_top$rank[[i]]),
      invitro_rank = as.integer(invitro_top$rank[[j]]),
      invivo_seed = as.integer(invivo_top$seed[[i]]),
      invitro_seed = as.integer(invitro_top$seed[[j]]),
      invivo_objective = as.numeric(invivo_top$objective[[i]]),
      invitro_objective = as.numeric(invitro_top$objective[[j]]),
      stringsAsFactors = FALSE
    )
    for (p in PARAMETERS) row[[paste0("ratio_", p)]] <- ratio_vals[[match(p, PARAMETERS)]]
    for (p in PARAMETERS) row[[paste0("log10_ratio_", p)]] <- log_ratio_vals[[match(p, PARAMETERS)]]
    pair_rows[[idx]] <- row
    idx <- idx + 1L
  }
}
pairs <- do.call(rbind, pair_rows)

matrix_cols <- paste0("log10_ratio_", PARAMETERS)
x <- as.matrix(pairs[, matrix_cols, drop = FALSE])
colnames(x) <- PARAMETERS

set.seed(umap_seed)
embedding <- uwot::umap(
  x,
  n_neighbors = min(15L, nrow(x) - 1L),
  min_dist = 0.1,
  metric = "euclidean",
  n_components = 2L,
  n_threads = 1L,
  ret_model = FALSE,
  verbose = FALSE
)
coords <- cbind(
  pairs[, c("pair_id", "invivo_rank", "invitro_rank", "invivo_seed", "invitro_seed", "invivo_objective", "invitro_objective")],
  data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], stringsAsFactors = FALSE)
)

plot_df <- coords
plot_df$invitro_rank_factor <- factor(plot_df$invitro_rank, levels = seq_len(top_n))

shape_values <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4)
names(shape_values) <- as.character(seq_along(shape_values))

subtitle <- paste(
  paste0(nrow(pairs), " pairings; 14 dimensions from log10(ratio_vivo_to_vitro):"),
  paste(strwrap(paste(PARAMETERS, collapse = ", "), width = 115), collapse = "\n"),
  sep = "\n"
)

p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(
    aes(color = invivo_rank, shape = invitro_rank_factor),
    size = 3.8,
    stroke = 1.05
  ) +
  geom_text_repel(
    aes(label = pair_id),
    color = "grey15",
    size = 3.0,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.alpha = 0.35,
    box.padding = 0.25,
    point.padding = 0.15,
    seed = umap_seed
  ) +
  scale_color_gradient(
    low = "#2b6cb0",
    high = "#d33f3f",
    limits = c(1, top_n),
    breaks = unique(round(seq(1, top_n, length.out = 5), 1)),
    name = "In Vivo Rank"
  ) +
  scale_shape_manual(values = shape_values[as.character(seq_len(top_n))], name = "In Vitro Rank") +
  guides(
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)
  ) +
  labs(
    title = "UMAP of Separate-Fit In Vivo/In Vitro Ratio Warm Starts",
    subtitle = subtitle,
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 22, face = "plain"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.7),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ratio_path <- file.path(out_dir, "cross_paired_top10_ratio_matrix.tsv")
coords_path <- file.path(out_dir, "cross_paired_top10_umap_coords.tsv")
png_path <- file.path(out_dir, "joint_soft_coupling_ratio_umap_500seed.png")
pdf_path <- file.path(out_dir, "joint_soft_coupling_ratio_umap_500seed.pdf")

write_tsv(pairs, ratio_path)
write_tsv(coords, coords_path)
ggsave(png_path, p, width = 16, height = 10, dpi = 220, bg = "white")
ggsave(pdf_path, p, width = 16, height = 10, bg = "white")

message("Wrote ratio matrix: ", ratio_path)
message("Wrote UMAP coordinates: ", coords_path)
message("Wrote UMAP PNG: ", png_path)
message("Wrote UMAP PDF: ", pdf_path)
message("Top in vivo seeds: ", paste(invivo_top$seed, collapse = ","))
message("Top in vitro seeds: ", paste(invitro_top$seed, collapse = ","))
