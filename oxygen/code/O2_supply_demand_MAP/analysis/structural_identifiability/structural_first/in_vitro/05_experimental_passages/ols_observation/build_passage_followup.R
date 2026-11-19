# List all documented transfers on the 114 OLS passages' physical ancestry chains.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script_path, mustWork = TRUE))
repo <- normalizePath(file.path(here, "../../../.."), mustWork = TRUE)
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else here
path <- file.path(
  dirname(here), "20261117_culture_actual_schedule_step01", "passage_crosswalk.tsv"
)
snapshot <- file.path(
  dirname(repo), "soft_couping_org/data/InVitroData",
  "cloneid_passaging_sum159_snapshot_20260731.tsv"
)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
crosswalk <- read_tsv(path)
raw <- read_tsv(snapshot)
stopifnot(!anyDuplicated(crosswalk$passage_id), !anyDuplicated(raw$id))
ols_index <- which(is.finite(crosswalk$recorded_ols_growth_per_day))
ancestry <- function(index) {
  chain <- integer()
  while (!is.na(index)) {
    stopifnot(!index %in% chain)
    chain <- c(chain, index)
    index <- match(crosswalk$physical_parent_passage_id[[index]], crosswalk$passage_id)
  }
  chain
}
needed <- sort(unique(unlist(lapply(ols_index, ancestry))))
rows <- crosswalk[needed, ]
seeds <- raw[match(rows$passage_id, raw$id), ]
donors <- raw[match(rows$direct_parent_record_id, raw$id), ]
stopifnot(!anyNA(seeds$id), !anyNA(donors$id), all(seeds$event == "seeding"),
          all(donors$event == "harvest"),
          all(seeds$passaged_from_id1 == donors$id),
          all(seeds$date == rows$seeding_time), all(donors$date == rows$direct_parent_record_time))
seed_at <- as.POSIXct(seeds$date, tz = "UTC")
donor_at <- as.POSIXct(donors$date, tz = "UTC")
days <- as.numeric(difftime(seed_at, donor_at, units = "days"))
stopifnot(all(days >= 0), all(abs(days - rows$recorded_parent_to_seed_days) < 1e-10))
followup <- data.frame(
  cohort = rows$cohort,
  child_passage_id = seeds$id,
  child_is_ols_passage = needed %in% ols_index,
  child_is_analysis_root = is.na(rows$physical_parent_passage_id),
  donor_record_id = donors$id,
  donor_passage_id_in_crosswalk = rows$physical_parent_passage_id,
  donor_timestamp_as_recorded = donors$date,
  child_seed_timestamp_as_recorded = seeds$date,
  recorded_gap_days = sprintf("%.9f", days),
  recorded_gap_over_7_days_for_review = days > 7,
  donor_comment_as_recorded = donors$comment,
  seed_comment_as_recorded = seeds$comment
)
stopifnot(nrow(followup) == 118L, sum(followup$child_is_ols_passage) == 114L,
          sum(followup$child_is_analysis_root) == 2L,
          sum(followup$recorded_gap_over_7_days_for_review) == 13L)
followup <- followup[order(-days, followup$child_passage_id), ]
write.table(followup, file.path(out, "passage_gap_followup.tsv"), sep = "\t",
            row.names = FALSE, quote = TRUE, na = "")
cat(sprintf(
  "%d candidate transfer records, including %d upstream passages and %d recorded gaps over 7 days.\n",
  nrow(followup), sum(!followup$child_is_ols_passage),
  sum(followup$recorded_gap_over_7_days_for_review)
))
