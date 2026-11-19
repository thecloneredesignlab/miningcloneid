# Run from the repository root with Rscript --vanilla. Uses archived trajectories only.
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
out <- "oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step01"
bridge <- "oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01"
contract <- "oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03"
anchor <- "oxygen/results/structural_identifiability/20261109_seed144_wp0/source/seed144"
snapshot_path <- "../soft_couping_org/data/InVitroData/cloneid_passaging_sum159_snapshot_20260731.tsv"
passages <- read_tsv("data/InVitroData_LTEE/invitro_passage_observations.tsv")
snapshot <- read_tsv(snapshot_path)
metadata <- read.csv("oxygen/data/metadata.csv")
raw_manifest <- read_tsv(file.path(contract, "growth_raw_count_manifest.tsv"))
kary <- read_tsv("data/InVitroData_LTEE/invitro_kary_cells.tsv")
schedule <- read_tsv(file.path(contract, "experimental_schedule.tsv"))
replay <- read_tsv(file.path(bridge, "verification/R01/selection.tsv"))
counts <- read_tsv(file.path(anchor, "invitro_daily_counts.tsv"))
summary <- read_tsv(file.path(anchor, "invitro_lineage_summary.tsv"))
id <- passages$passage_id
stopifnot(!anyDuplicated(id), !anyDuplicated(snapshot$id), !anyDuplicated(raw_manifest$record_id),
          !anyDuplicated(metadata$id), setequal(summary$passage_id, id))

# Load function definitions without initialising a model or running source-file entry points.
prod <- new.env()
for (part in c("lineage_utils", "objective_utils", "io_utils", "summary_utils", "lineage_simulation_utils")) {
  file <- file.path(bridge, "inputs/production", paste0("o2_supply_demand_map_invitro_", part, ".R"))
  for (expr in parse(file)) {
    if (is.call(expr) && as.character(expr[[1L]]) == "<-" && is.symbol(expr[[2L]]) &&
        is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
      eval(expr, prod)
    }
  }
}
fit_dir <- file.path(bridge, "inputs/data")
fit <- readRDS(file.path(fit_dir, "fit_data.Rds"))
fit <- prod$ivt_attach_g0g1_flow_data(fit, fit_dir, file.path(fit_dir, "g0g1_ploidy_density_grid.csv"))
stopifnot(setequal(names(fit), id))
cohorts <- c("2N", "4N")
adapters <- lapply(cohorts, function(cohort) {
  jobs <- readRDS(file.path(fit_dir, paste0("jobs_", cohort, ".Rds")))
  prod$ivt_build_job_table_adapter(jobs, fit, cohort)$segments
})
segments <- unlist(adapters, recursive = FALSE)
key <- function(cohort, segment) paste(cohort, segment)
names(segments) <- vapply(segments, function(s) key(s$cohort, s$segment_id), "")
segment_by_id <- unlist(lapply(segments, function(s) setNames(rep(key(s$cohort, s$segment_id),
                                                               length(s$data_ids)), s$data_ids)),
                        use.names = TRUE)
names(segment_by_id) <- unlist(lapply(segments, function(s) s$data_ids), use.names = FALSE)
stopifnot(!anyDuplicated(names(segments)), !anyDuplicated(names(segment_by_id)),
          setequal(names(segment_by_id), id), identical(names(segments), key(schedule$cohort, schedule$segment_id)))
stopifnot(!anyDuplicated(key(replay$cohort, replay$segment_id)),
          setequal(names(segments), key(replay$cohort, replay$segment_id)))
replay <- replay[match(names(segments), key(replay$cohort, replay$segment_id)), ]
counts_by_segment <- split(counts, key(counts$cohort, counts$segment_id))
segment_for_passage <- segments[unname(segment_by_id[id])]
summary <- summary[match(id, summary$passage_id), ]
parent_for_passage <- vapply(segment_for_passage, function(s) s$parent_segment_id, "")
stopifnot(all(passages$oxygen_pct == vapply(segment_for_passage, function(s) s$oxygen_pct, 0)),
          identical(unname(parent_for_passage), summary$parent_segment_id),
          identical(unname(vapply(segments, function(s) s$duration_days, 0)), as.numeric(schedule$nominal_duration_days)),
          all(vapply(segments, function(s) identical(as.numeric(s$obs_days_local),
              as.numeric(strsplit(schedule$candidate_days[match(key(s$cohort, s$segment_id),
                  names(segments))], ";", fixed = TRUE)[[1L]])), FALSE)))

# Reapply original selectors to archived daily live totals. This checks indices and
# growth, while state-composition inheritance is verified from the production source.
initial_total <- readRDS(file.path(fit_dir, "fit_result.rds"))$cfg$init_total_size
selection <- vector("list", length(segments))
rates <- numeric(length(segments))
for (i in seq_along(segments)) {
  s <- segments[[i]]
  target <- s$final_cells
  if ((!is.finite(target) || target <= 0) && i < length(segments) &&
      segments[[i + 1L]]$cohort == s$cohort) target <- segments[[i + 1L]]$initial_cells
  requested <- if (is.finite(s$initial_cells) && s$initial_cells > 0) s$initial_cells else initial_total
  daily <- counts_by_segment[[names(segments)[[i]]]]
  daily <- daily[order(daily$day), ]
  sim <- list(Ntot_live_obs = daily$live_cells, live_state_obs = matrix(daily$live_cells, ncol = 1))
  stopifnot(all(daily$day == s$obs_days_local), all(daily$dead_total_cells[1L] == 0))
  picked <- prod$ivt_select_segment_observation(sim, requested, 1, target, s$obs_days_local)
  stopifnot(isTRUE(all.equal(picked$target_live_cells, schedule$target_live_cells[[i]],
                             tolerance = 1e-12)),
            picked$selected_index == replay$index[[i]], picked$selected_day == replay$day[[i]])
  if (is.na(s$parent_segment_id)) {
    stopifnot(is.na(replay$parent_day[[i]]))
  } else {
    parent <- selection[[match(key(s$cohort, s$parent_segment_id), names(segments))]]
    boundary <- prod$ivt_extract_passage_end_state(parent$sim, requested, 1, parent$target, parent$days)
    stopifnot(boundary$passage_executed, boundary$selected_index == replay$parent_index[[i]],
              boundary$selected_day == replay$parent_day[[i]])
  }
  rates[[i]] <- prod$ivt_log_growth_rate(daily$live_cells[[1L]], picked$selected_live_cells, picked$selected_day)
  selection[[i]] <- list(sim = sim, target = target, days = s$obs_days_local)
}
i <- match(unname(segment_by_id[id]), names(segments))
stopifnot(max(abs(rates[i] - summary$predicted_growth_rate)) < 1e-10,
          all(summary$selected_day == replay$day[i]),
          all(summary$selected_day == schedule$selected_day[i]))

# Resolve the nearest included seeding ancestor on the database primary-parent chain.
seeds <- snapshot[match(id, snapshot$id), ]
stopifnot(all(seeds$id == id), all(seeds$event == "seeding"), all(is.na(seeds$passaged_from_id2) | !nzchar(seeds$passaged_from_id2)))
parent_id <- function(seed) {
  current <- seed
  visited <- character()
  repeat {
    parent <- snapshot$passaged_from_id1[match(current, snapshot$id)]
    if (is.na(parent) || !nzchar(parent)) return(NA_character_)
    stopifnot(!parent %in% visited, parent != seed)
    if (parent %in% id) return(parent)
    visited <- c(visited, parent)
    current <- parent
  }
}
parents <- vapply(id, parent_id, "")
has_parent <- !is.na(parents)
stopifnot(sum(has_parent) == 129L, all(has_parent == !is.na(parent_for_passage)),
          all(segment_by_id[parents[has_parent]] == key(passages$cohort[has_parent],
                                                       parent_for_passage[has_parent])))
direct_parents <- snapshot[match(seeds$passaged_from_id1, snapshot$id), ]
stopifnot(all(!is.na(direct_parents$id)), all(direct_parents$event == "harvest"))
gap <- as.numeric(difftime(as.POSIXct(seeds$date, tz = "UTC"),
                           as.POSIXct(direct_parents$date, tz = "UTC"), units = "days"))
observed_media <- snapshot$media[match(metadata$id, snapshot$id)]
media_labels <- unique(data.frame(media = observed_media, oxygen = metadata$oxygen_pct))
stopifnot(!anyNA(media_labels), !anyDuplicated(media_labels$media))
oxygen_from_media <- media_labels$oxygen[match(seeds$media, media_labels$media)]
stopifnot(all(is.na(oxygen_from_media) | oxygen_from_media == passages$oxygen_pct))

# Use the exact seeding-plus-direct-children grouping and the original date-based OLS.
manifest_counts <- snapshot[match(raw_manifest$record_id, snapshot$id), ]
stopifnot(all(raw_manifest$record_id == manifest_counts$id),
          all(raw_manifest$date == as.Date(manifest_counts$date)),
          all(raw_manifest$corrected_count == manifest_counts$correctedCount),
          all(raw_manifest$physical_parent == manifest_counts$passaged_from_id1))
observed_growth <- which(is.finite(passages$observed_growth))
stopifnot(setequal(id[observed_growth], raw_manifest$passage_id))
metadata_in_manifest <- raw_manifest$source == "oxygen/data/metadata.csv"
m <- metadata[match(raw_manifest$record_id[metadata_in_manifest], metadata$id), ]
stopifnot(all(m$id == raw_manifest$record_id[metadata_in_manifest]),
          all(m$date == raw_manifest$date[metadata_in_manifest]),
          all(m$correctedCount == raw_manifest$corrected_count[metadata_in_manifest]))
count_groups <- split(raw_manifest, raw_manifest$passage_id)
ols <- endpoint <- rep(NA_real_, length(id))
ols_days <- ols_count_ids <- ols_count_timestamps <- rep(NA_character_, length(id))
for (j in observed_growth) {
  p <- id[[j]]
  raw <- snapshot[(snapshot$id == p | snapshot$passaged_from_id1 %in% p) &
                  is.finite(snapshot$correctedCount), ]
  raw <- raw[order(raw$date), ]
  day <- as.numeric(as.Date(raw$date) - min(as.Date(raw$date)))
  log_count <- log(pmax(1, raw$correctedCount))
  stopifnot(all(raw$id %in% count_groups[[p]]$record_id),
            all(count_groups[[p]]$record_id %in% raw$id),
            max(day) == passages$passage_duration[[j]],
            abs(raw$correctedCount[[1L]] - passages$initial_cells[[j]]) < 0.5,
            abs(tail(raw$correctedCount, 1L) - passages$final_cells[[j]]) < 0.5)
  ols[[j]] <- unname(coef(lm(log_count ~ day))[[2L]])
  endpoint[[j]] <- (tail(log_count, 1L) - log_count[[1L]]) / max(day)
  ols_days[[j]] <- paste(day, collapse = ";")
  ols_count_ids[[j]] <- paste(raw$id, collapse = ";")
  ols_count_timestamps[[j]] <- paste(raw$date, collapse = ";")
}
ols_error <- abs(ols - passages$observed_growth)
stopifnot(max(ols_error, na.rm = TRUE) < 1e-12)
for (j in seq_along(id)) {
  obs <- prod$ivt_observed_passage_summary(fit[[id[[j]]]])
  observed_cells <- kary$observed_kary_N[kary$passage_id == id[[j]]]
  stopifnot(identical(sort(obs$observed_kary), sort(as.numeric(observed_cells))),
            obs$observed_n_kary == passages$n_kary[[j]], obs$observed_n_flow == passages$n_flow[[j]])
  if (obs$observed_n_flow > 0L) {
    stopifnot(obs$observed_flow_sample_name == passages$flow_sample_name[[j]],
              all(diff(obs$observed_flow$ploidy) > 0),
              abs(sum(obs$observed_flow$density) * median(diff(obs$observed_flow$ploidy)) - 1) < 1e-10)
  }
}

harvests <- snapshot[snapshot$event == "harvest" & snapshot$passaged_from_id1 %in% id, ]
harvests <- harvests[order(harvests$date), ]
harvest_groups <- split(harvests, harvests$passaged_from_id1)[id]
stopifnot(all(vapply(harvest_groups, nrow, integer(1)) > 0L))
all_harvest_dates <- vapply(harvest_groups, function(d) paste(d$date, collapse = ";"), "")
last_harvest <- vapply(harvest_groups, function(d) tail(d$date, 1L), "")
seed_to_last_harvest_day <- as.numeric(as.Date(last_harvest) - as.Date(seeds$date))
fallback <- vapply(segments, function(s) all(vapply(s$observed, function(o) is.na(o$passage_duration), FALSE)), FALSE)
crosswalk <- data.frame(
  passage_id = id, cohort = passages$cohort, physical_parent_passage_id = unname(parents),
  direct_parent_record_id = direct_parents$id, direct_parent_record_time = direct_parents$date,
  seeding_time = seeds$date, recorded_parent_to_seed_days = gap, linked_harvest_times = unname(all_harvest_dates),
  seed_to_last_recorded_harvest_calendar_days = seed_to_last_harvest_day,
  media_id = seeds$media, assigned_oxygen_pct = passages$oxygen_pct,
  oxygen_pct_from_count_metadata_media_join = oxygen_from_media,
  seed_corrected_count = seeds$correctedCount, seed_raw_cell_count = seeds$cellCount,
  figure_cumulative_day = passages$cumulative_end_day,
  production_segment_id = summary$segment_id, production_parent_segment_id = summary$parent_segment_id,
  production_nominal_duration_days = schedule$nominal_duration_days[i],
  production_uses_fallback_duration = unname(fallback[i]),
  seed144_selected_local_day = replay$day[i], seed144_parent_reseed_local_day = replay$parent_day[i],
  seed144_endpoint_growth_per_day = rates[i], ols_count_calendar_days = ols_days,
  ols_count_record_ids = ols_count_ids, ols_count_timestamps = ols_count_timestamps,
  recorded_ols_growth_per_day = passages$observed_growth, reconstructed_ols_per_day = ols,
  observed_first_last_count_rate_per_day = endpoint, karyotype_cell_count = passages$n_kary,
  flow_sample_name = passages$flow_sample_name, flow_processed_grid_points = passages$n_flow)
write.table(crosswalk, file.path(out, "passage_crosswalk.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)

cat("Passages", nrow(crosswalk), "by cohort", table(crosswalk$cohort), "linked ancestry edges", sum(has_parent), "\n")
cat("Production segments", lengths(adapters), "tree edges", sum(!is.na(replay$parent_day)),
    "fallback durations", sum(fallback),
    "endpoint vs nominal changes", sum(replay$day != schedule$nominal_duration_days),
    "no positive target", sum(is.na(schedule$target_live_cells)),
    "growth placements with selected vs OLS duration change",
    sum(replay$day[i[observed_growth]] != passages$passage_duration[observed_growth]), "\n")
cat("Parent-to-child model propagation and parent observation day discrepancies",
    sum(replay$parent_day[i[has_parent]] != replay$day[match(segment_by_id[parents[has_parent]], names(segments))]), "\n")
cat("Oxygen labels corroborated by media/count join", sum(!is.na(oxygen_from_media)),
    "missing media mappings", sum(is.na(oxygen_from_media)), "\n")
cat("No fitted duration", sum(is.na(passages$passage_duration)), "with dated harvest range",
    range(seed_to_last_harvest_day[is.na(passages$passage_duration)]),
    "days. Parent-to-seed intervals >=30 days", sum(gap[has_parent] >= 30), "\n")
cat("OLS groups", length(observed_growth), "native count rows", nrow(raw_manifest),
    "count records per passage", range(vapply(count_groups, nrow, integer(1))), "max OLS error", max(ols_error, na.rm = TRUE),
    "OLS vs same-count endpoint abs difference >1e-8",
    sum(abs(ols - endpoint) > 1e-8, na.rm = TRUE), "max difference",
    max(abs(ols - endpoint), na.rm = TRUE), "\n")
cat("Anchor endpoint-growth max replay discrepancy", max(abs(rates[i] - summary$predicted_growth_rate)), "\n")
cat("Karyotype samples", sum(passages$n_kary > 0), "cell counts", sum(passages$n_kary),
    "processed flow curves", sum(passages$n_flow > 0), "grids", unique(passages$n_flow[passages$n_flow > 0]), "\n")
