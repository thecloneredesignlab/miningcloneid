# Shared, side-effect-free aggregation for the 500 fitted in-vivo endpoints.
suppressPackageStartupMessages(library(data.table))

supp4_4_metrics <- c("ploidy", "effective_o2_pct", "tumor_size_mm3")
supp4_4_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
supp4_4_quantiles <- function(x) {
  supp4_4_assert(length(x) > 0L && all(is.finite(x)), "Non-finite/empty seed distribution.")
  unname(stats::quantile(x, c(0.25, 0.5, 0.75), type = 7L))
}

supp4_4_aggregate <- function(tumor, expected_seeds = seq_len(500L)) {
  tumor <- as.data.table(copy(tumor))
  keys <- c("seed_number", "harvest", "cohort", "day")
  supp4_4_assert(all(c(keys, supp4_4_metrics) %in% names(tumor)), "Missing trajectory columns.")
  supp4_4_assert(!anyDuplicated(tumor[, ..keys]), "Duplicate seed/tumor/day rows.")
  supp4_4_assert(setequal(unique(tumor$seed_number), expected_seeds), "Seed set is incomplete.")
  supp4_4_assert(all(is.finite(as.matrix(tumor[, ..supp4_4_metrics]))), "Non-finite trajectory values.")
  coverage <- tumor[, .(n_seeds = uniqueN(seed_number)), by = .(harvest, cohort, day)]
  supp4_4_assert(all(coverage$n_seeds == length(expected_seeds)), "Missing seed at a tumor/time point.")
  cohort <- tumor[, c(
    list(n_tumors = uniqueN(harvest)),
    lapply(.SD, mean)
  ), by = .(seed_number, cohort, day), .SDcols = supp4_4_metrics]
  setorder(cohort, seed_number, cohort, day)
  cohort[, segment := rleid(n_tumors), by = .(seed_number, cohort)]
  long <- melt(cohort, id.vars = c("seed_number", "cohort", "day", "n_tumors", "segment"),
               measure.vars = supp4_4_metrics, variable.name = "metric", value.name = "value")
  summary <- long[, {
    q <- supp4_4_quantiles(value)
    list(n_seeds = uniqueN(seed_number), n_tumors = unique(n_tumors),
         segment = unique(segment), q25 = q[1L], median = q[2L], q75 = q[3L])
  }, by = .(cohort, day, metric)]
  supp4_4_assert(all(summary$n_seeds == length(expected_seeds)), "Summary does not contain all seeds.")
  supp4_4_assert(!anyDuplicated(summary[, .(cohort, day, metric)]), "Inconsistent cohort membership across seeds.")
  supp4_4_assert(all(summary$q25 <= summary$median & summary$median <= summary$q75), "Invalid quantile ordering.")
  setorder(summary, cohort, metric, day)
  list(cohort = cohort, summary = summary)
}
