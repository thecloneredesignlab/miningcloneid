# Recheck the actual observations' binary64 kernel and pruning masks in R 4.4.2.
bridge <- "oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01"
out <- "oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step02"
stored <- readRDS(file.path(bridge, "verification/R01/result.rds"))
p <- stored$params
masks <- read.delim(file.path(out, "anchor_384/probability_masks.tsv"), colClasses = "character")
grid <- 22:154
stopifnot(getRversion() == "4.4.2", .Machine$double.min.exp == -1022, .Machine$double.digits == 53)
for (i in seq_len(nrow(masks))) {
  o <- as.numeric(masks$oxygen_pct[[i]])
  n <- as.integer(masks$chr_count[[i]])
  hill <- p$O2_crit^p$n_O / (p$O2_crit^p$n_O + o^p$n_O)
  death <- p$mu_hp * hill * (n / 44)^p$gamma_mu
  mis <- p$p_mis_base + p$p_misseg * death / (p$k_o_mis + death)
  keep <- paste0(as.integer(dbinom(0:n, n, mis) >= 1e-8), collapse = "")
  stopifnot(identical(keep, masks$kept_mask[[i]]))
}

prod <- new.env()
objective <- file.path(bridge, "inputs/production/o2_supply_demand_map_invitro_objective_utils.R")
for (expr in parse(objective)) {
  if (is.call(expr) && identical(expr[[2L]], as.name("ivt_smooth_kary_distribution"))) eval(expr, prod)
}
for (width in seq(.014, .015, length.out=65)) {
  kernel <- outer(grid, grid, function(obs_n, true_n) stats::dnorm(obs_n, true_n, width))
  normalized <- sweep(kernel, 2, colSums(kernel), "/", check.margin = FALSE)
  stopifnot(identical(unname(normalized), diag(length(grid))))
  stopifnot(all(vapply(stored$pred$kary, function(prob) {
    identical(prod$ivt_smooth_kary_distribution(grid, prob, width),
              prod$ivt_smooth_kary_distribution(grid, prob, p$sigma_kary))
  }, FALSE)))
}
cat("R 4.4.2 binary64\n", nrow(masks), "oxygen-by-mother masks all match\n",
    "65 widths produced identity kernels and bitwise-equal predictions for 12 archived distributions\n",
    sep = "")

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  # Export the independent native daily state reference for the outward-ball replay review.
  ids <- stored$pred$selection
  items <- lapply(seq_len(nrow(ids)), function(i) {
    key <- paste(ids$cohort[[i]], ids$segment_id[[i]])
    a <- stored$pred$states[[key]]
    data.frame(cohort=ids$cohort[[i]], segment=ids$segment_id[[i]], day=seq_len(nrow(a))-1L, a)
  })
  options(digits = 17)
  write.table(do.call(rbind, items), args[[1L]], sep="\t", row.names=FALSE, quote=FALSE)
}
