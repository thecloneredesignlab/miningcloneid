args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript plot_joint.R OUTPUT_DIR DESIGN_DIR ARCHIVED_FIT_DIR")
directory <- normalizePath(args[1], mustWork = TRUE)
design <- normalizePath(args[2], mustWork = TRUE)
fit <- normalizePath(args[3], mustWork = TRUE)
if (file.exists(file.path(directory, "SHA256SUMS"))) stop("Choose a fresh output directory for the plot")
read_tsv <- function(path) read.delim(path, check.names = FALSE, na.strings = "")
canonical <- function(value) sub("SUM159_NLS", "SUM-159_NLS", value, fixed = TRUE)

prediction <- subset(read_tsv(file.path(fit, "readout_predictions.tsv")), mode == "primary")
prediction$event_id <- canonical(prediction$event_id)
mask <- read_tsv(file.path(design, "event_mask.tsv"))
specimens <- read_tsv(file.path(design, "specimen_decisions.tsv"))
seeds <- subset(specimens, decision == "conditioned_seed_record")
seeds$date <- as.Date(seeds$timestamp_as_recorded)
identifiers <- match(paste(prediction$readout, prediction$event_id),
                     paste(mask$modality, mask$event))
stopifnot(length(identifiers) == 426, !anyNA(identifiers))
prediction$date <- as.Date(mask$timestamp[identifiers])
prediction$path <- paste(prediction$cohort, prediction$lineage)
prediction$branch <- grepl("_A7K_", prediction$event_id, fixed = TRUE)
cohorts <- c("2N", "4N")
lineages <- c("C", "O1", "O2")
obs <- "#152D46"
pred <- "#AC4E36"
latent <- "#087B70"

pdf(file.path(directory, "joint_trajectories.pdf"), width = 12, height = 8.7,
    title = "Joint experimental culture observations and archived predictions")
draw_calendar <- function(readout, title, axis_label, factor, show_live = FALSE) {
  d <- prediction[prediction$readout == readout, ]
  par(mfrow = c(2, 3), mar = c(3.6, 4.1, 3.1, 1), oma = c(2.4, 0, 3.2, 0))
  for (cohort in cohorts) {
    for (line in lineages) {
      b <- d[d$cohort == cohort & d$lineage == line, ]
      if (!nrow(b)) {
        plot.new()
        title(paste(cohort, line))
        text(.5, .5, "No classified observations")
        if (cohort == "2N" && line == "C") {
          legend("topleft", c("Observed", "Conditional joint fit"), col = c(obs, pred),
                 pch = c(16, 1), bty = "n", cex = .84)
        }
        next
      }
      ys <- c(b$observed, b$predicted, if (show_live) b$predicted_live_objects else NULL) * factor
      transfer <- seeds$date[seeds$cohort == cohort & seeds$lineage == line]
      calendar_range <- range(c(b$date, transfer))
      plot(b$date, b$observed * factor, ylim = range(c(0, ys)), type = "n",
           xlim = calendar_range, main = paste(cohort, line), ylab = axis_label, xlab = "",
           xaxt = "n", bty = "l", cex.lab = .87)
      marks <- pretty(calendar_range, n = 4)
      axis.Date(1, at = marks, format = "%b %d", cex.axis = .78)
      rug(as.numeric(transfer), ticksize = .018, col = "#7B8790")
      arrows(b$date, b$observed * factor, b$date, b$predicted * factor,
             length = 0, col = adjustcolor("grey45", .45), lwd = .8)
      if (show_live) points(b$date, b$predicted_live_objects * factor, pch = 3, col = latent, cex = .7)
      size <- if (show_live) rep(.83, nrow(b)) else pmax(.54, pmin(1.16, sqrt(b$denominator / 7000)))
      points(b$date, b$observed * factor, pch = ifelse(b$branch, 18, 16), col = obs, cex = size)
      points(b$date, b$predicted * factor, pch = ifelse(b$branch, 5, 1), col = pred, cex = size)
      if (cohort == "2N" && line == "C") {
        labs <- c("Observed", "Joint fit", if (show_live) "Model live component")
        colors <- c(obs, pred, if (show_live) latent)
        legend("topleft", legend = labs, col = colors, pch = c(16, 1, if (show_live) 3),
               bty = "n", cex = .82)
      }
    }
  }
  mtext(title, side = 3, outer = TRUE, line = 1.15, font = 2, cex = 1.3)
  mtext(if (show_live) "Dates show recorded harvests. Diamonds mark A7K. Short ticks mark seeded passages."
        else "Diamonds mark A7K. Point area reflects eligible-object denominators. Short ticks mark seeded passages.",
        side = 1, outer = TRUE, line = .65, cex = .84)
}

draw_calendar("corrected_count", "Corrected visible-object counts and model live component",
              "Objects (millions)", 1e-6, TRUE)
draw_calendar("death_fraction", "Measured death fractions and conditional joint predictions",
              "Classified dead (%)", 100)

kary <- subset(read_tsv(file.path(fit, "karyotype_distributions.tsv")), mode == "primary")
kary$event_id <- canonical(kary$event_id)
chromosome_events <- prediction[prediction$readout == "karyotype", ]
stopifnot(nrow(kary) == 798, sum(kary$observed_cells) == 104)
par(mfrow = c(2, 3), mar = c(4, 4.25, 3.2, 1), oma = c(2.4, 0, 3.2, 0))
for (i in seq_len(nrow(chromosome_events))) {
  e <- chromosome_events[i, ]
  d <- kary[kary$event_id == e$event_id, ]
  n <- sum(d$observed_cells)
  y <- d$observed_cells/n
  label <- paste(e$cohort, e$lineage, if (e$branch) "A7K" else "A17")
  plot(d$autosomes, y, type = "n", ylim = range(0, y, d$fitted_probability),
       xlab = "Autosome count", ylab = "Probability", main = label, bty = "l")
  segments(d$autosomes, 0, d$autosomes, y, col = adjustcolor(obs, .78), lwd = 2)
  lines(d$autosomes, d$fitted_probability, col = pred, lwd = 1.7)
  mtext(paste(format(e$date, "%d %b %Y"), "|", n, "cells"), side = 3, line = .15, cex = .74)
  if (i == 1) legend("topright", c("Observed", "Joint fit"), col = c(obs, pred), lwd = 2, bty = "n", cex = .84)
}
mtext("Full observed and predicted chromosome distributions", side = 3, outer = TRUE,
      line = 1.15, font = 2, cex = 1.3)
mtext("Each panel retains all 133 autosome bins. Counts and predictions use the matched harvest event.",
      side = 1, outer = TRUE, line = .65, cex = .84)
dev.off()
