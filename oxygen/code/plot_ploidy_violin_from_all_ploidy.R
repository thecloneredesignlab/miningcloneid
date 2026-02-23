#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0) return(out)
  for (a in argv) {
    if (!startsWith(a, "--")) next
    kv <- sub("^--", "", a)
    if (grepl("=", kv, fixed = TRUE)) {
      k <- sub("=.*$", "", kv)
      v <- sub("^[^=]*=", "", kv)
      out[[k]] <- v
    } else {
      out[[kv]] <- TRUE
    }
  }
  out
}

dose_to_colors <- function(dose, low_col, high_col) {
  d <- suppressWarnings(as.numeric(dose))
  t <- rep(0.5, length(d))
  ok <- is.finite(d)
  if (sum(ok) >= 2) {
    dmin <- min(d[ok], na.rm = TRUE)
    dmax <- max(d[ok], na.rm = TRUE)
    if (isTRUE(dmax > dmin)) {
      t[ok] <- (d[ok] - dmin) / (dmax - dmin)
    }
  } else if (sum(ok) == 1) {
    t[ok] <- 1.0
  }
  t <- pmax(0, pmin(1, t))
  rgb_mat <- grDevices::colorRamp(c(low_col, high_col), space = "Lab")(t)
  grDevices::rgb(rgb_mat[, 1], rgb_mat[, 2], rgb_mat[, 3], maxColorValue = 255)
}

build_violin_plot <- function(dat, sample_info, title, subtitle) {
  sample_levels <- sample_info$sample
  dat$sample <- factor(dat$sample, levels = sample_levels)
  dat$ploidy_group <- sample_info$ploidy_group[match(as.character(dat$sample), sample_info$sample)]
  dat$ploidy_group_plot <- ifelse(dat$ploidy_group %in% c("2N", "4N"), dat$ploidy_group, "Other")
  dat$ploidy_group_plot <- factor(dat$ploidy_group_plot, levels = c("2N", "4N", "Other"))

  violin_fill <- setNames(sample_info$fill_color, sample_info$sample)
  line_colors <- c("2N" = "#2E8B57", "4N" = "#C0392B", "Other" = "#4D4D4D")

  ggplot(dat, aes(x = sample, y = ploidy, fill = sample, color = ploidy_group_plot)) +
    geom_violin(linewidth = 0.30, scale = "width", trim = TRUE) +
    geom_boxplot(aes(group = sample), width = 0.10, outlier.size = 0.2, alpha = 0.65, fill = "white", color = "black", linewidth = 0.25) +
    stat_summary(aes(group = sample), fun = median, geom = "point", shape = 21, size = 1.3, fill = "white", color = "black", stroke = 0.2) +
    scale_fill_manual(values = violin_fill, guide = "none") +
    scale_color_manual(values = line_colors, breaks = c("2N", "4N"), drop = FALSE, name = "Init ploidy") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample",
      y = "Ploidy"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

main <- function() {
  argv <- parse_args(commandArgs(trailingOnly = TRUE))

  in_tsv <- argv$in_tsv %||%
    "/Users/4482173/Documents/GitHub/miningcloneid/data/InVivoData_Gemcitabine/all_ploidy.tsv"
  out_pdf <- argv$out_pdf %||%
    "/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/ploidy_violin_by_sample.pdf"
  out_pdf_dose0 <- argv$out_pdf_dose0 %||%
    sub("\\.pdf$", "_dose0.pdf", out_pdf)
  if (identical(out_pdf_dose0, out_pdf)) {
    out_pdf_dose0 <- paste0(out_pdf, "_dose0.pdf")
  }

  if (!file.exists(in_tsv)) {
    stop("Input file not found: ", in_tsv)
  }

  dat <- read.delim(in_tsv, check.names = FALSE, stringsAsFactors = FALSE)
  req <- c("file", "ploidy")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  dat <- dat[is.finite(dat$ploidy), , drop = FALSE]
  if (nrow(dat) == 0) stop("No valid ploidy values after filtering.")

  dat$sample <- sub("\\.sps\\.cbs$", "", dat$file)

  sample_info <- data.frame(sample = unique(dat$sample), stringsAsFactors = FALSE)
  parts <- strsplit(sample_info$sample, "-", fixed = TRUE)
  sample_info$ploidy_group <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  sample_info$dose <- suppressWarnings(as.numeric(vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))))

  sample_info$ploidy_order <- ifelse(sample_info$ploidy_group == "2N", 1L, ifelse(sample_info$ploidy_group == "4N", 2L, 99L))
  sample_info <- sample_info[order(sample_info$ploidy_order, sample_info$dose, sample_info$sample), , drop = FALSE]

  sample_info$fill_color <- "#BDBDBD"
  idx_2n <- sample_info$ploidy_group == "2N"
  idx_4n <- sample_info$ploidy_group == "4N"
  if (any(idx_2n)) {
    sample_info$fill_color[idx_2n] <- dose_to_colors(sample_info$dose[idx_2n], low_col = "#D9F0D3", high_col = "#006D2C")
  }
  if (any(idx_4n)) {
    sample_info$fill_color[idx_4n] <- dose_to_colors(sample_info$dose[idx_4n], low_col = "#FDD0C7", high_col = "#A50F15")
  }

  p <- build_violin_plot(
    dat = dat,
    sample_info = sample_info,
    title = "Ploidy Distribution by Sample",
    subtitle = "Ordered by ploidy (2N -> 4N) and dose (low -> high); darker color means higher dose"
  )

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_pdf, p, width = 12, height = 6.5, device = grDevices::cairo_pdf)

  dose0_samples <- sample_info[is.finite(sample_info$dose) & sample_info$dose == 0, , drop = FALSE]
  if (nrow(dose0_samples) > 0) {
    dat_dose0 <- dat[as.character(dat$sample) %in% dose0_samples$sample, , drop = FALSE]
    p_dose0 <- build_violin_plot(
      dat = dat_dose0,
      sample_info = dose0_samples,
      title = "Ploidy Distribution by Sample (Dose = 0)",
      subtitle = "Dose-zero samples only; ordered by ploidy (2N -> 4N)"
    )
    dir.create(dirname(out_pdf_dose0), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_pdf_dose0, p_dose0, width = 10, height = 6, device = grDevices::cairo_pdf)
    message("Done. Outputs written to: ", normalizePath(out_pdf), " ; ", normalizePath(out_pdf_dose0))
  } else {
    message("Done. Output written to: ", normalizePath(out_pdf), " (no dose=0 samples found for extra plot)")
  }
}

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
}

if (sys.nframe() == 0) {
  main()
}
