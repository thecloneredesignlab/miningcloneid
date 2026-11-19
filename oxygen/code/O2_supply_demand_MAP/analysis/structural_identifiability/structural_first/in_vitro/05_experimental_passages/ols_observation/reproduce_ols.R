# Reproduce the experimental slopes directly from the recorded CloneID counts.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script_path, mustWork = TRUE))
repo <- normalizePath(file.path(here, "../../../.."), mustWork = TRUE)
results <- dirname(here)
step1 <- file.path(results, "20261118_culture_experimental_OLS_step01")
snapshot <- file.path(
  dirname(repo), "soft_couping_org/data/InVitroData",
  "cloneid_passaging_sum159_snapshot_20260731.tsv"
)
args <- commandArgs(TRUE)
out <- if (length(args)) args[[1L]] else here
dir.create(out, showWarnings = FALSE, recursive = TRUE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
mapped <- read_tsv(file.path(step1, "count_event_mapping.tsv"))
raw <- read_tsv(snapshot)
obs <- read_tsv(file.path(repo, "data/InVitroData_LTEE/invitro_passage_observations.tsv"))
obs <- obs[is.finite(obs$observed_growth), ]
stopifnot(nrow(mapped) == 335L, !anyDuplicated(mapped$record_id),
          nrow(obs) == 114L, !anyDuplicated(obs$passage_id), !anyDuplicated(raw$id),
          setequal(mapped$passage_id, obs$passage_id))

counts <- raw[match(mapped$record_id, raw$id), ]
stopifnot(identical(counts$id, mapped$record_id),
          all(counts$correctedCount == mapped$db_correctedCount),
          all(counts$date == mapped$count_timestamp_as_recorded),
          all(is.finite(counts$correctedCount)))
counts$passage_id <- mapped$passage_id
counts$calendar_date <- as.Date(counts$date, format = "%Y-%m-%d")
stopifnot(!anyNA(counts$calendar_date))

slopes <- lapply(split(counts, counts$passage_id), function(rows) {
  rows <- rows[order(rows$calendar_date, rows$id), ]
  passage_id <- rows$passage_id[[1L]]
  reference <- obs[match(passage_id, obs$passage_id), ]
  days <- as.integer(rows$calendar_date - min(rows$calendar_date))
  stopifnot(length(days) >= 2L, !anyDuplicated(days),
            rows$id[[1L]] == passage_id, rows$event[[1L]] == "seeding",
            all(rows$event[-1L] == "harvest"),
            all(rows$passaged_from_id1[-1L] == passage_id),
            max(days) == reference$passage_duration,
            all(days == mapped$ols_calendar_day[match(rows$id, mapped$record_id)]))
  fit <- lm(log(pmax(1, rows$correctedCount)) ~ days)
  slope <- unname(coef(fit)[["days"]])
  data.frame(
    passage_id = passage_id,
    cohort = reference$cohort,
    n_counts = nrow(rows),
    calendar_days = paste(days, collapse = ";"),
    recorded_ols_per_day = reference$observed_growth,
    reproduced_ols_per_day = slope,
    abs_difference_per_day = abs(slope - reference$observed_growth)
  )
})
slopes <- do.call(rbind, slopes)
tolerance <- 1e-12
stopifnot(nrow(slopes) == 114L, sum(slopes$n_counts) == 335L,
          all(slopes$abs_difference_per_day <= tolerance))
numeric_fields <- c("recorded_ols_per_day", "reproduced_ols_per_day", "abs_difference_per_day")
slopes_out <- slopes
slopes_out[numeric_fields] <- lapply(slopes[numeric_fields], function(x) sprintf("%.17g", x))
write.table(slopes_out, file.path(out, "experimental_ols_reproduction.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf(
  "%d counts, %d reproduced slopes, maximum absolute error %.17g/day (tolerance %.1g).\n",
  sum(slopes$n_counts), nrow(slopes), max(slopes$abs_difference_per_day), tolerance
))
