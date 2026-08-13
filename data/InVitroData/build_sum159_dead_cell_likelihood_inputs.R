#!/usr/bin/env Rscript

# Build model-aware SUM-159 dead-cell likelihood inputs using only local files.
# This script intentionally does not load cloneid or connect to the database.

options(stringsAsFactors = FALSE, digits = 15)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve this script's path.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
data_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(data_dir, "..", ".."), mustWork = TRUE)

snapshot_date <- "20260731"
dead_path <- file.path(data_dir, "sum159_all_db_cloneid_dead_cell_frequency.tsv")
db_source_path <- file.path(
  data_dir,
  paste0("cloneid_passaging_dead_cell_sources_and_parents_", snapshot_date, ".tsv")
)
db_sum159_path <- file.path(
  data_dir,
  paste0("cloneid_passaging_sum159_snapshot_", snapshot_date, ".tsv")
)
db_metadata_path <- file.path(
  data_dir,
  paste0("cloneid_passaging_snapshot_metadata_", snapshot_date, ".tsv")
)
fit_data_path <- file.path(repo_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "fit_data.Rds")
jobs_2n_path <- file.path(repo_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "jobs_2N.Rds")
jobs_4n_path <- file.path(repo_root, "oxygen", "ploidyOxygen", "data", "fit_objects", "jobs_4N.Rds")

required_paths <- c(
  dead_path, db_source_path, db_sum159_path, db_metadata_path,
  fit_data_path, jobs_2n_path, jobs_4n_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required local inputs: ", paste(missing_paths, collapse = ", "))
}

read_tsv <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character()
  )
}

write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

canonicalize_id <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  sub("^SUM159_", "SUM-159_", x)
}

parse_formal_passage_index <- function(x) {
  x <- as.character(x)
  hit <- grepl("_(O1|O2)_A[0-9]+_seed$", x)
  out <- rep(NA_integer_, length(x))
  out[hit] <- as.integer(sub("^.*_A([0-9]+)_seed$", "\\1", x[hit]))
  out
}

dead <- read_tsv(dead_path)
db_source <- read_tsv(db_source_path)
db_sum159 <- read_tsv(db_sum159_path)
db_metadata <- read_tsv(db_metadata_path)
fit_data <- readRDS(fit_data_path)
jobs_2n <- readRDS(jobs_2n_path)
jobs_4n <- readRDS(jobs_4n_path)

required_dead_columns <- c(
  "db_source_id", "source_id", "context_key", "ploidy", "condition_replicate",
  "event", "biological_passage", "ltee_comment_membership",
  "classification_status", "n_bundles_success", "n_bundles_failed",
  "n_bundles_unavailable", "failure_reasons", "cell_state", "count",
  "eligible_denominator", "frequency", "feature_modes", "prediction_shard"
)
missing_dead_columns <- setdiff(required_dead_columns, names(dead))
if (length(missing_dead_columns) > 0L) {
  stop("Dead-cell TSV is missing columns: ", paste(missing_dead_columns, collapse = ", "))
}
if (nrow(dead) != 240L || anyDuplicated(dead$db_source_id)) {
  stop("Expected 240 unique dead-cell database source IDs.")
}
if (!identical(as.character(dead$db_source_id), as.character(db_source$db_source_id))) {
  stop("Local DB source snapshot no longer aligns row-for-row with the dead-cell TSV.")
}
if (!all(dead$db_source_id %in% db_sum159$id)) {
  stop("Local SUM-159 DB snapshot does not contain every dead-cell source ID.")
}
if (!identical(as.character(dead$event), as.character(db_source$db_event))) {
  stop("Event mismatch between dead-cell TSV and local DB snapshot.")
}

dead_count <- suppressWarnings(as.numeric(dead$count))
eligible_denominator <- suppressWarnings(as.numeric(dead$eligible_denominator))
observed_dead_fraction <- suppressWarnings(as.numeric(dead$frequency))
if (
  any(!is.finite(dead_count)) || any(dead_count < 0) ||
  any(!is.finite(eligible_denominator)) || any(eligible_denominator <= 0) ||
  any(dead_count > eligible_denominator) ||
  any(!is.finite(observed_dead_fraction)) ||
  max(abs(observed_dead_fraction - dead_count / eligible_denominator)) > 1e-8
) {
  stop("Dead-cell count, denominator, or frequency validation failed.")
}

fit_ids <- names(fit_data)
fit_ids_o1o2 <- fit_ids[grepl("_(O1|O2)_", fit_ids)]
if (length(fit_ids_o1o2) != 94L) {
  stop("Expected 94 O1/O2 fit-data identifiers; found ", length(fit_ids_o1o2), ".")
}

flatten_jobs <- function(jobs, cohort) {
  rows <- lapply(seq_len(nrow(jobs)), function(i) {
    ids <- as.character(jobs$data_ids[[i]])
    data.frame(
      model_passage_id = ids,
      cohort_from_jobs = cohort,
      source_sim_key = rep(as.character(jobs$sim_key[[i]]), length(ids)),
      source_parent_key = rep(as.character(jobs$parent_key[[i]]), length(ids)),
      source_depth = rep(as.integer(jobs$depth[[i]]), length(ids)),
      model_oxygen_pct = rep(as.numeric(jobs$oxygen[[i]]), length(ids)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

job_map <- rbind(flatten_jobs(jobs_2n, "2N"), flatten_jobs(jobs_4n, "4N"))
if (anyDuplicated(job_map$model_passage_id)) {
  stop("A fit-data identifier occurs in more than one jobs row.")
}

choose_model_passage_id <- function(i) {
  candidates <- if (identical(as.character(dead$event[[i]]), "seeding")) {
    canonicalize_id(dead$db_source_id[[i]])
  } else {
    canonicalize_id(c(db_source$passaged_from_id1[[i]], db_source$passaged_from_id2[[i]]))
  }
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[candidates %in% fit_ids_o1o2]
  if (length(hits) > 1L) {
    stop("Ambiguous model passage mapping for ", dead$db_source_id[[i]], ": ", paste(hits, collapse = ", "))
  }
  if (length(hits) == 0L) NA_character_ else hits[[1]]
}

model_passage_id <- vapply(seq_len(nrow(dead)), choose_model_passage_id, character(1))
model_match <- match(model_passage_id, fit_ids)
job_match <- match(model_passage_id, job_map$model_passage_id)

model_passage_duration <- rep(NA_real_, nrow(dead))
model_initial_cells <- rep(NA_real_, nrow(dead))
model_final_cells <- rep(NA_real_, nrow(dead))
mapped_idx <- which(!is.na(model_match))
model_passage_duration[mapped_idx] <- vapply(
  fit_data[model_match[mapped_idx]],
  function(x) suppressWarnings(as.numeric(x$passage_duration)),
  numeric(1)
)
model_initial_cells[mapped_idx] <- vapply(
  fit_data[model_match[mapped_idx]],
  function(x) suppressWarnings(as.numeric(x$initial_cells)),
  numeric(1)
)
model_final_cells[mapped_idx] <- vapply(
  fit_data[model_match[mapped_idx]],
  function(x) suppressWarnings(as.numeric(x$final_cells)),
  numeric(1)
)

lineage_passage_index <- parse_formal_passage_index(model_passage_id)
is_formal_main_chain <- !is.na(lineage_passage_index) &
  is.finite(model_passage_duration) & model_passage_duration > 0
is_a7k_side_branch <- !is.na(model_passage_id) & grepl("_A7K_seed$", model_passage_id)
model_mapping_status <- ifelse(
  is_formal_main_chain,
  "formal_main_chain",
  ifelse(is_a7k_side_branch, "a7k_side_branch_not_currently_simulated", "not_in_current_fit_objects")
)

parent1_is_seed <- as.character(db_source$parent1_event) == "seeding"
parent2_is_seed <- as.character(db_source$parent2_event) == "seeding"
is_seeding <- as.character(dead$event) == "seeding"
is_harvest <- as.character(dead$event) == "harvest"
seed_db_id <- ifelse(
  is_seeding,
  as.character(dead$db_source_id),
  ifelse(parent1_is_seed, as.character(db_source$passaged_from_id1),
         ifelse(parent2_is_seed, as.character(db_source$passaged_from_id2), NA_character_))
)
seed_datetime <- ifelse(
  is_seeding,
  as.character(db_source$db_event_datetime),
  ifelse(parent1_is_seed, as.character(db_source$parent1_datetime),
         ifelse(parent2_is_seed, as.character(db_source$parent2_datetime), NA_character_))
)
elapsed_days_since_seed <- ifelse(
  is_seeding,
  0,
  ifelse(parent1_is_seed, suppressWarnings(as.numeric(db_source$elapsed_days_from_parent1)),
         ifelse(parent2_is_seed, suppressWarnings(as.numeric(db_source$elapsed_days_from_parent2)), NA_real_))
)
if (any(is_harvest & !is.finite(elapsed_days_since_seed))) {
  stop(
    "Harvest rows without a valid seeding-parent elapsed time: ",
    paste(dead$db_source_id[is_harvest & !is.finite(elapsed_days_since_seed)], collapse = ", ")
  )
}
if (any(is_harvest & elapsed_days_since_seed < 0)) {
  stop("At least one harvest precedes its seeding parent.")
}

endpoint_offset_days <- ifelse(
  is_harvest & is_formal_main_chain,
  elapsed_days_since_seed - model_passage_duration,
  NA_real_
)
abs_endpoint_offset_days <- abs(endpoint_offset_days)
endpoint_selection_rank <- rep(NA_integer_, nrow(dead))
selected_for_endpoint <- rep(FALSE, nrow(dead))

formal_harvest_idx <- which(is_harvest & is_formal_main_chain)
formal_by_passage <- split(formal_harvest_idx, model_passage_id[formal_harvest_idx])
for (pid in names(formal_by_passage)) {
  idx <- formal_by_passage[[pid]]
  ord <- order(
    abs_endpoint_offset_days[idx],
    as.character(db_source$db_event_datetime[idx]),
    as.character(dead$db_source_id[idx])
  )
  endpoint_selection_rank[idx[ord]] <- seq_along(ord)
  selected_for_endpoint[idx[ord[[1]]]] <- TRUE
}

landmark_harvest_idx <- which(is_harvest & is_a7k_side_branch)
landmark_by_passage <- split(landmark_harvest_idx, model_passage_id[landmark_harvest_idx])
if (length(landmark_by_passage) != 4L || any(lengths(landmark_by_passage) != 1L)) {
  stop("Expected one A7K harvest observation for each of four O1/O2 cohort-scenarios.")
}
endpoint_selection_rank[landmark_harvest_idx] <- 1L
selected_for_endpoint[landmark_harvest_idx] <- TRUE

within_model_time_window <- is_harvest & is_formal_main_chain &
  is.finite(elapsed_days_since_seed) & elapsed_days_since_seed >= 0 &
  elapsed_days_since_seed <= model_passage_duration
include_in_current_endpoint_likelihood <- selected_for_endpoint & is_formal_main_chain
include_in_future_timecourse_likelihood <- within_model_time_window

scenario_id <- paste(dead$ploidy, dead$condition_replicate, sep = "-")
model_segment_id <- ifelse(
  is_formal_main_chain,
  paste0(scenario_id, "-A", lineage_passage_index),
  ifelse(is_a7k_side_branch, paste0(scenario_id, "-A7K"), NA_character_)
)
observation_role <- ifelse(
  is_seeding,
  "seeding_qc_only",
  ifelse(
    !is.na(model_passage_id) & is_a7k_side_branch,
    "a7k_side_branch_endpoint_candidate",
    ifelse(
      is.na(model_passage_id),
      "out_of_current_model",
      ifelse(selected_for_endpoint, "formal_fixed_endpoint", "formal_intra_passage")
    )
  )
)
exclusion_reason <- ifelse(
  include_in_current_endpoint_likelihood,
  "",
  ifelse(
    is_seeding,
    "seeding_event_current_model_resets_dead_stock_to_zero",
    ifelse(
      is.na(model_passage_id),
      "source_parent_not_present_in_current_fit_objects",
      ifelse(
        is_a7k_side_branch,
        "a7k_side_branch_not_simulated_by_current_independent_scenario_adapter",
        "not_closest_harvest_to_fixed_endpoint"
      )
    )
  )
)

metadata_lookup <- setNames(as.character(db_metadata$value), as.character(db_metadata$key))
db_snapshot_server_time <- unname(metadata_lookup[["database_server_time"]])
db_snapshot_session_time_zone <- unname(metadata_lookup[["database_session_time_zone"]])

enriched <- data.frame(
  observation_id = canonicalize_id(dead$source_id),
  db_source_id = as.character(dead$db_source_id),
  source_id = as.character(dead$source_id),
  context_key = as.character(dead$context_key),
  cohort = as.character(dead$ploidy),
  ploidy = as.character(dead$ploidy),
  lineage_id = as.character(dead$condition_replicate),
  scenario_id = scenario_id,
  condition_replicate = as.character(dead$condition_replicate),
  lineage_group = "deprived",
  model_passage_id = model_passage_id,
  model_segment_id = model_segment_id,
  lineage_passage_index = lineage_passage_index,
  landmark_after_passage_index = ifelse(is_a7k_side_branch, 7L, NA_integer_),
  model_source_depth = ifelse(!is.na(job_match), job_map$source_depth[job_match], NA_integer_),
  model_oxygen_pct = ifelse(!is.na(job_match), job_map$model_oxygen_pct[job_match], NA_real_),
  model_mapping_status = model_mapping_status,
  observation_role = observation_role,
  event = as.character(dead$event),
  db_event = as.character(db_source$db_event),
  source_biological_passage = suppressWarnings(as.integer(dead$biological_passage)),
  db_biological_passage = suppressWarnings(as.integer(db_source$db_biological_passage)),
  db_event_datetime = as.character(db_source$db_event_datetime),
  seed_db_id = seed_db_id,
  seed_datetime = seed_datetime,
  elapsed_days_since_seed = elapsed_days_since_seed,
  model_endpoint_day = model_passage_duration,
  endpoint_offset_days = endpoint_offset_days,
  abs_endpoint_offset_days = abs_endpoint_offset_days,
  endpoint_selection_rank = endpoint_selection_rank,
  selected_for_endpoint = selected_for_endpoint,
  likelihood_observation_day = ifelse(
    include_in_current_endpoint_likelihood,
    model_passage_duration,
    NA_real_
  ),
  within_model_time_window = within_model_time_window,
  include_in_current_endpoint_likelihood = include_in_current_endpoint_likelihood,
  include_in_future_timecourse_likelihood = include_in_future_timecourse_likelihood,
  exclusion_reason = exclusion_reason,
  cell_state = as.character(dead$cell_state),
  dead_count = dead_count,
  eligible_denominator = eligible_denominator,
  observed_dead_fraction = observed_dead_fraction,
  classification_status = as.character(dead$classification_status),
  classification_partial = as.character(dead$classification_status) != "classified_full_promoted",
  n_bundles_success = suppressWarnings(as.integer(dead$n_bundles_success)),
  n_bundles_failed = suppressWarnings(as.integer(dead$n_bundles_failed)),
  n_bundles_unavailable = suppressWarnings(as.integer(dead$n_bundles_unavailable)),
  failure_reasons = as.character(dead$failure_reasons),
  ltee_comment_membership = as.logical(dead$ltee_comment_membership),
  feature_modes = as.character(dead$feature_modes),
  prediction_shard = as.character(dead$prediction_shard),
  model_initial_cells = model_initial_cells,
  model_final_cells = model_final_cells,
  db_snapshot_file = basename(db_source_path),
  db_snapshot_server_time = db_snapshot_server_time,
  db_snapshot_session_time_zone = db_snapshot_session_time_zone,
  stringsAsFactors = FALSE
)

harvest_timecourse <- enriched[
  enriched$event == "harvest" & !is.na(enriched$model_passage_id),
  ,
  drop = FALSE
]
endpoint_candidates <- harvest_timecourse[harvest_timecourse$selected_for_endpoint, , drop = FALSE]
endpoint_ready <- endpoint_candidates[
  endpoint_candidates$include_in_current_endpoint_likelihood,
  ,
  drop = FALSE
]

sort_rows <- function(x) {
  x[order(
    factor(x$cohort, levels = c("2N", "4N")),
    factor(x$lineage_id, levels = c("O1", "O2")),
    ifelse(is.na(x$lineage_passage_index), Inf, x$lineage_passage_index),
    x$db_event_datetime,
    x$db_source_id
  ), , drop = FALSE]
}
harvest_timecourse <- sort_rows(harvest_timecourse)
endpoint_candidates <- sort_rows(endpoint_candidates)
endpoint_ready <- sort_rows(endpoint_ready)

if (nrow(harvest_timecourse) != 199L) {
  stop("Expected 199 current-fit mapped harvest observations; found ", nrow(harvest_timecourse), ".")
}
if (nrow(endpoint_candidates) != 94L || anyDuplicated(endpoint_candidates$model_passage_id)) {
  stop("Expected 94 unique formal/A7K endpoint candidates.")
}
if (nrow(endpoint_ready) != 90L || anyDuplicated(endpoint_ready$model_passage_id)) {
  stop("Expected 90 unique current-model formal endpoint likelihood rows.")
}
if (!setequal(endpoint_ready$model_passage_id, fit_ids_o1o2[!grepl("_A7K_seed$", fit_ids_o1o2)])) {
  stop("Endpoint-ready rows do not exactly cover the 90 formal O1/O2 fit-data passages.")
}
if (any(endpoint_ready$event != "harvest") || any(!endpoint_ready$selected_for_endpoint)) {
  stop("Endpoint-ready output contains a non-harvest or non-selected observation.")
}
if (max(endpoint_ready$abs_endpoint_offset_days) > 0.5) {
  stop("A selected formal endpoint is more than 0.5 days from the fixed model endpoint.")
}

output_paths <- c(
  db_enriched = file.path(data_dir, paste0("sum159_dead_cell_db_enriched_", snapshot_date, ".tsv")),
  harvest_timecourse = file.path(data_dir, paste0("sum159_dead_cell_harvest_timecourse_", snapshot_date, ".tsv")),
  endpoint_candidates = file.path(data_dir, paste0("sum159_dead_cell_endpoint_candidates_", snapshot_date, ".tsv")),
  endpoint_ready = file.path(data_dir, paste0("sum159_dead_cell_endpoint_likelihood_ready_", snapshot_date, ".tsv"))
)
write_tsv(enriched, output_paths[["db_enriched"]])
write_tsv(harvest_timecourse, output_paths[["harvest_timecourse"]])
write_tsv(endpoint_candidates, output_paths[["endpoint_candidates"]])
write_tsv(endpoint_ready, output_paths[["endpoint_ready"]])

ready_counts <- as.data.frame.matrix(table(endpoint_ready$cohort, endpoint_ready$lineage_id))
summary_metrics <- c(
  source_dead_cell_rows = nrow(dead),
  db_enriched_rows = nrow(enriched),
  harvest_rows_total = sum(enriched$event == "harvest"),
  harvest_rows_mapped_to_current_fit_objects = nrow(harvest_timecourse),
  endpoint_candidates_including_a7k = nrow(endpoint_candidates),
  endpoint_ready_formal_main_chain = nrow(endpoint_ready),
  endpoint_candidates_a7k_excluded_current_model = sum(endpoint_candidates$model_mapping_status == "a7k_side_branch_not_currently_simulated"),
  out_of_current_model_rows = sum(enriched$model_mapping_status == "not_in_current_fit_objects"),
  seeding_rows_qc_only = sum(enriched$event == "seeding"),
  endpoint_ready_partial_classification_rows = sum(endpoint_ready$classification_partial),
  endpoint_ready_max_abs_time_offset_days = max(endpoint_ready$abs_endpoint_offset_days),
  endpoint_ready_median_abs_time_offset_days = median(endpoint_ready$abs_endpoint_offset_days),
  endpoint_ready_2N_O1 = ready_counts["2N", "O1"],
  endpoint_ready_2N_O2 = ready_counts["2N", "O2"],
  endpoint_ready_4N_O1 = ready_counts["4N", "O1"],
  endpoint_ready_4N_O2 = ready_counts["4N", "O2"]
)
summary_df <- data.frame(
  metric = names(summary_metrics),
  value = as.character(unname(summary_metrics)),
  stringsAsFactors = FALSE
)
summary_path <- file.path(data_dir, paste0("sum159_dead_cell_likelihood_build_summary_", snapshot_date, ".tsv"))
write_tsv(summary_df, summary_path)

manifest_files <- c(
  dead_path, db_source_path, db_sum159_path, db_metadata_path,
  fit_data_path, jobs_2n_path, jobs_4n_path, script_path,
  unname(output_paths), summary_path
)
manifest <- data.frame(
  file = basename(manifest_files),
  absolute_path = normalizePath(manifest_files, mustWork = TRUE),
  md5 = unname(tools::md5sum(manifest_files)),
  bytes = as.numeric(file.info(manifest_files)$size),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(data_dir, paste0("sum159_dead_cell_likelihood_manifest_", snapshot_date, ".tsv"))
write_tsv(manifest, manifest_path)

cat("DB_ENRICHED=", output_paths[["db_enriched"]], "\n", sep = "")
cat("HARVEST_TIMECOURSE=", output_paths[["harvest_timecourse"]], "\n", sep = "")
cat("ENDPOINT_CANDIDATES=", output_paths[["endpoint_candidates"]], "\n", sep = "")
cat("ENDPOINT_READY=", output_paths[["endpoint_ready"]], "\n", sep = "")
cat("SUMMARY=", summary_path, "\n", sep = "")
cat("MANIFEST=", manifest_path, "\n", sep = "")
