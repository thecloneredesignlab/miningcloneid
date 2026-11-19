# Generate the calendar-day sampling/weight table and verify it by observed-count replay.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script_path, mustWork = TRUE))
results <- dirname(here)
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else here
source(file.path(here, "observation_operator.R"))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
counts <- read_tsv(file.path(
  results, "20261118_culture_experimental_OLS_step01", "count_event_mapping.tsv"
))
fit <- read_tsv(file.path(
  results, "20261109_seed144_wp0/source/seed144", "fit_summary.tsv"
))
slopes <- read_tsv(file.path(out, "experimental_ols_reproduction.tsv"))
mapping <- build_ols_calendar_map(counts, as.numeric(fit$value[fit$metric == "dt"]))

# Replay with recorded visible counts to check the linear combination independently of lm().
replay <- mapping[c("cohort", "model_segment_id", "calendar_day")]
replay$visible_count <- counts$db_correctedCount[match(mapping$record_id, counts$record_id)]
replayed <- apply_ols_calendar_map(mapping, replay)
stopifnot(identical(replayed$passage_id, slopes$passage_id))
max_error <- max(abs(replayed$predicted_ols_per_day - slopes$reproduced_ols_per_day))

weights_ok <- vapply(split(mapping, mapping$passage_id), function(rows) {
  w <- rows$ols_weight_numerator
  denominator <- unique(rows$ols_weight_denominator)
  length(denominator) == 1L && sum(w) == 0 && sum(w * rows$calendar_day) == denominator
}, logical(1))
stopifnot(nrow(mapping) == 335L, nrow(replayed) == 114L, all(weights_ok),
          sum(mapping$observed_after_seed144_selected_day) == 33L, max_error <= 1e-12)
write.table(mapping, file.path(out, "ols_calendar_map.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat(sprintf(
  "Calendar map: %d samples, %d slopes, maximum replay error vs lm() %.17g/day.\n",
  nrow(mapping), nrow(replayed), max_error
))
