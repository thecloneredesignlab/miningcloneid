# Run from the StructuralAndPracticalIdentifiability root with Rscript --vanilla.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else getwd())
syn <- file.path(root, "oxygen/results/synthetic_identifiability")
struct <- file.path(root, "oxygen/results/structural_identifiability")
out <- file.path(syn, "20261119_experimental_karyotype_alignment_01")
input_repo <- file.path(dirname(root), "ploidyOxygen")
raw_dir <- file.path(dirname(root), "soft_couping_org/data/InVitroData")
hist <- file.path(struct, "20261118_culture_experimental_inheritance_01/hist")
counts_dir <- file.path(syn, "20261118_experimental_counts_practical_01")
read_tsv <- function(path) {
  read.delim(path, stringsAsFactors = FALSE, colClasses = "character")
}
write_tsv <- function(data, name) {
  write.table(data, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
}
find_one <- function(data, key, value) {
  row <- data[!is.na(data[[key]]) & data[[key]] == value, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  row
}
as_flag <- function(value) tolower(value) == "true"
collapse <- function(value) paste(value, collapse = "|")
record_date <- function(x) substr(x, 1L, 10L)
age_day <- function(record, seed) as.integer(as.Date(record_date(record)) -
                                              as.Date(record_date(seed)))
stamp <- function(x) format(x, format = "%Y-%m-%d %H:%M:%S")
main <- function() {
  commit <- "5c4fbe5bc0885ba9d2e3a095d7abe74c727980c5"
  bundle_file <- "data/workflow_inputs/workflow_inputs.Rds"
  bundle_path <- tempfile(fileext = ".rds")
  on.exit(unlink(bundle_path), add = TRUE)
  git_arg <- function(command, path) {
    c("-C", shQuote(input_repo), command, paste0(commit, ":", path))
  }
  stopifnot(system2("git", git_arg("show", bundle_file), stdout = bundle_path) == 0L,
            identical(system2("git", c("hash-object", shQuote(bundle_path)),
                              stdout = TRUE), "0c65880f7cc27d1c1bd405a5096af73583cb5371"))
  bundle <- readRDS(bundle_path)
  fit_path <- file.path(root, "oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds")
  stopifnot(identical(system2("git", c("hash-object", shQuote(fit_path)), stdout = TRUE),
                      system2("git", git_arg("rev-parse", "data/fit_objects/fit_data.Rds"),
                              stdout = TRUE)))
  fit <- readRDS(fit_path)
  fit <- lapply(Filter(function(x) length(x$kary) > 0L, fit), `[[`, "kary")

  data_dir <- file.path(root, "data/InVitroData_LTEE")
  cells <- read_tsv(file.path(data_dir, "invitro_kary_cells.tsv"))
  archived <- read_tsv(file.path(struct,
    "20261109_seed144_wp0/source/seed144/invitro_observed_kary.tsv"))
  contract <- read_tsv(file.path(syn,
    "20261113_culture_contract_steps01_03/observed_karyotype_cells.tsv"))
  passage_obs <- read_tsv(file.path(data_dir, "invitro_passage_observations.tsv"))
  count <- read_tsv(file.path(counts_dir, "count_inputs.tsv"))
  graph <- read_tsv(file.path(hist, "passage_design.tsv"))
  proxy <- read_tsv(file.path(hist, "dp_proxy_histograms.tsv"))
  restarts <- read_tsv(file.path(hist, "restart_policy.tsv"))
  raw <- read_tsv(file.path(raw_dir, "cloneid_passaging_sum159_snapshot_20260731.tsv"))
  death <- read_tsv(file.path(raw_dir, "sum159_dead_cell_harvest_timecourse_20260731.tsv"))
  stopifnot(nrow(cells) == 220L, length(bundle$cloneIds) == 12L, length(fit) == 12L,
            !anyDuplicated(cells[c("passage_id", "cell_index")]),
            nrow(count) == 335L, !anyDuplicated(count$record_id),
            nrow(graph) == 116L, !anyDuplicated(graph$passage_id),
            sum(as_flag(count$fitted_count)) == 221L)
  common <- c("passage_id", "cell_index", "observed_kary_N")
  names(contract) <- common
  canonical <- function(data, fields) {
    data <- data[order(data$passage_id, as.integer(data$cell_index)), fields]
    row.names(data) <- NULL
    data
  }
  archived_fields <- setdiff(names(archived), c("passage_index", "segment_id"))
  stopifnot(identical(canonical(cells, common), canonical(contract, common)),
            identical(canonical(cells, archived_fields),
                      canonical(archived, archived_fields)))

  raw_forms <- lapply(bundle$cloneIds, function(origin) {
    mat <- do.call(rbind, bundle$karyotypes$karyotype[bundle$karyotypes$id == origin])
    model_n <- rowSums(mat[, !grepl("^9999:", colnames(mat)), drop = FALSE])
    autosomes <- grepl("^([1-9]|1[0-9]|2[0-2]):", colnames(mat))
    stopifnot(sum(autosomes) == 22L)
    fit_key <- names(fit)[vapply(fit, identical, logical(1L), y = model_n)]
    stopifnot(length(fit_key) == 1L)
    data.frame(fit_passage_id = fit_key, assay_origin_id = origin,
               cell_index = seq_along(model_n),
               observed_N_archived = model_n,
               observed_N_22_autosomes = rowSums(mat[, autosomes, drop = FALSE]),
               observed_N_locus_999 = mat[, grepl("^999:", colnames(mat))])
  })
  raw_forms <- do.call(rbind, raw_forms)
  stopifnot(!anyDuplicated(raw_forms[c("fit_passage_id", "cell_index")]),
            nrow(raw_forms) == 220L,
            all(raw_forms$observed_N_22_autosomes + raw_forms$observed_N_locus_999 ==
                raw_forms$observed_N_archived))

  fitted <- count[as_flag(count$fitted_count), , drop = FALSE]
  summary <- lapply(unique(cells$passage_id), function(key) {
    cc <- cells[cells$passage_id == key, , drop = FALSE]
    cc <- cc[order(as.integer(cc$cell_index)), , drop = FALSE]
    x <- raw_forms[raw_forms$fit_passage_id == key, , drop = FALSE]
    origin <- unique(x$assay_origin_id)
    stopifnot(length(origin) == 1L, identical(as.integer(cc$cell_index), x$cell_index),
              identical(as.numeric(cc$observed_kary_N), x$observed_N_archived))
    old <- find_one(passage_obs, "passage_id", key)
    stopifnot(as.integer(old$n_kary) == nrow(x), all(cc$oxygen_pct == old$oxygen_pct),
              all(cc$segment_id == old$segment_id))
    event <- find_one(raw, "id", origin)
    seed_id <- if (event$event == "harvest") event$passaged_from_id1 else origin
    seed <- find_one(raw, "id", seed_id)
    stopifnot(seed$event == "seeding")
    for (rec in list(event, seed)) {
      ref <- find_one(bundle$passaging, "id", rec$id)
      stopifnot(identical(stamp(ref$date), rec$date), ref$event == rec$event,
                as.character(ref$media) == rec$media,
                ref$passaged_from_id1 == rec$passaged_from_id1)
    }
    actual_o2 <- bundle$media$oxygen_pct[match(event$media, bundle$media$id)]
    stopifnot(!is.na(actual_o2),
              actual_o2 == bundle$media$oxygen_pct[match(seed$media, bundle$media$id)])
    day <- age_day(event$date, seed$date)
    in_graph <- seed_id %in% graph$passage_id
    in_fitted <- event$id %in% fitted$record_id
    if (in_graph) {
      stopifnot(actual_o2 == as.numeric(find_one(graph, "passage_id", seed_id)$oxygen_pct))
    }
    if (in_fitted) {
      joined <- find_one(fitted, "record_id", event$id)
      stopifnot(joined$passage_id == seed_id, as.integer(joined$calendar_day) == day,
                as.numeric(joined$oxygen_pct) == actual_o2)
    }
    proxy_source <- key %in% proxy$archival_donor_passage_id
    proxy_donor <- ""
    if (proxy_source) {
      h <- proxy[proxy$archival_donor_passage_id == key, , drop = FALSE]
      stopifnot(all(as.integer(h$total_cells) == nrow(x)),
                identical(sort(rep(as.numeric(h$chromosome_count), as.integer(h$cells))),
                          sort(x$observed_N_archived)))
      rr <- restarts[restarts$new_policy == "independent_dp_proxy" &
                     restarts$donor_passage_id_in_crosswalk == key, , drop = FALSE]
      stopifnot(nrow(rr) == 1L, rr$donor_record_id != origin)
      proxy_donor <- rr$donor_record_id
    }
    physical_counts <- fitted[fitted$passage_id == seed_id, , drop = FALSE]
    d <- death[death$db_source_id == origin, , drop = FALSE]
    stopifnot(nrow(d) <= 1L)
    status <- if (proxy_source) "proxy_only" else if (in_fitted) {
      "direct_harvest_anchor"
    } else if (in_graph) {
      "seeding_id_assay_day_unverified"
    } else {
      "outside_stage2_graph"
    }
    data.frame(
      cohort = old$cohort, archived_fit_passage_id = key,
      archived_model_segment = old$segment_id, archived_nominal_o2_pct = old$oxygen_pct,
      archived_model_cumulative_end_day = old$cumulative_end_day,
      assay_origin_id = origin, physical_passage_id = seed_id,
      source_event = event$event, source_event_timestamp = event$date,
      physical_seed_timestamp = seed$date, seed_donor_record = seed$passaged_from_id1,
      seed_record_corrected_count = seed$correctedCount,
      origin_record_corrected_count = event$correctedCount,
      nominal_o2_pct = actual_o2, record_calendar_day = day,
      harvest_anchored_assay_calendar_day = if (event$event == "harvest") day else NA,
      n_cells = nrow(x), archived_N_min = min(x$observed_N_archived),
      archived_N_max = max(x$observed_N_archived),
      autosomal_N_min = min(x$observed_N_22_autosomes),
      autosomal_N_max = max(x$observed_N_22_autosomes),
      locus_999_range = collapse(range(x$observed_N_locus_999)),
      physical_passage_in_stage2 = in_graph,
      assay_id_in_221 = in_fitted,
      physical_passage_fitted_count_days = collapse(physical_counts$calendar_day),
      physical_passage_fitted_count_ids = collapse(physical_counts$record_id),
      dp_archival_donor_record = proxy_donor, role = status,
      dead_n = if (nrow(d)) d$dead_count else NA,
      dead_eligible_n = if (nrow(d)) d$eligible_denominator else NA
    )
  })
  summary <- do.call(rbind, summary)
  expected <- c(direct_harvest_anchor = 40L, outside_stage2_graph = 100L,
                proxy_only = 40L, seeding_id_assay_day_unverified = 40L)
  totals <- tapply(summary$n_cells, summary$role, sum)
  stopifnot(identical(names(totals), names(expected)), all(totals == expected),
            sum(summary$archived_nominal_o2_pct != summary$nominal_o2_pct) == 4L)
  print(with(summary, table(role)))
  print(totals)
  print(with(summary, tapply(n_cells, nominal_o2_pct, sum)))
  death_overlap <- death[death$db_source_id %in% fitted$record_id, , drop = FALSE]
  stopifnot(nrow(death_overlap) == 195L, !anyDuplicated(death_overlap$db_source_id),
            sum(as_flag(death_overlap$include_in_future_timecourse_likelihood)) == 153L,
            all(!is.na(summary$dead_n[summary$role == "direct_harvest_anchor"])))
  cat("Matched morphology records:", nrow(death_overlap),
      "; marked in the legacy time window:",
      sum(as_flag(death_overlap$include_in_future_timecourse_likelihood)), "\n")
  write_tsv(summary, "sample_alignment.tsv")
  write_tsv(raw_forms, "cell_alignment.tsv")
}
main()
