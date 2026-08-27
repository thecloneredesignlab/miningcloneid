#!/usr/bin/env Rscript

script_dir <- local({
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[[1L]]))) else getwd()
})
source(file.path(script_dir, "util", "runtime", "process_runner.R"))

data_Figure1 <- function() {
  destination_root <- file.path(DATA_ROOT, "Figure1")
  frozen_source <- LTEE_DATA_ROOT
  frozen_files <- sort(list.files(
    frozen_source, full.names = TRUE, recursive = FALSE
  ))
  frozen_files <- frozen_files[file.info(frozen_files)$isdir %in% FALSE]
  expected <- c(
    "invitro_kary_cells.tsv",
    "invitro_lineage_timeline.tsv",
    "invitro_passage_observations.tsv",
    "invivo_burden_long.tsv",
    "invivo_harvest_catalog.tsv",
    "invivo_ploidy_cells.tsv"
  )
  if (!identical(basename(frozen_files), sort(expected))) {
    stop("Figure 1 frozen-input set differs from the approved six-file contract.")
  }

  soft_coupling_root <- normalizePath(
    file.path(INVITRO_RESULT_ROOT, "..", "..", ".."),
    mustWork = TRUE
  )
  population_source <- file.path(
    soft_coupling_root,
    "data", "InVitroData",
    "cloneid_passaging_sum159_snapshot_20260731.tsv"
  )
  flow_source <- file.path(
    INVITRO_RESULT_ROOT, "seed10", "invitro_observed_flow.tsv"
  )
  require_files(
    c(population_source, flow_source),
    "Figure 1 measurement source"
  )

  frozen_destinations <- file.path(destination_root, basename(frozen_files))
  frozen_copies <- Map(copy_input, frozen_files, frozen_destinations)
  raw_sources <- c(population_source, flow_source)
  raw_destinations <- file.path(
    destination_root, "source_raw", basename(raw_sources)
  )
  raw_copies <- Map(copy_input, raw_sources, raw_destinations)

  source_files <- c(frozen_files, raw_sources)
  copied <- c(frozen_copies, raw_copies)
  contract <- data.frame(
    role = c(
      rep("approved Figure 1 frozen input", length(frozen_files)),
      "raw population-size audit source",
      "observed flow-cytometry density source"
    ),
    source = normalizePath(source_files, mustWork = TRUE),
    local_file = unlist(copied, use.names = FALSE),
    source_md5 = unname(tools::md5sum(source_files)),
    local_md5 = unname(tools::md5sum(unlist(copied, use.names = FALSE))),
    stringsAsFactors = FALSE
  )
  write_data_contract("Figure1", contract)

  passage <- utils::read.delim(
    file.path(destination_root, "invitro_passage_observations.tsv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  population_snapshot <- utils::read.delim(
    raw_destinations[[1L]], check.names = FALSE, stringsAsFactors = FALSE
  )
  observed_flow <- utils::read.delim(
    raw_destinations[[2L]], check.names = FALSE, stringsAsFactors = FALSE
  )

  passage_ids <- as.character(passage$passage_id)
  if (anyDuplicated(passage_ids)) {
    stop("Figure 1 passage IDs must be unique before population-size auditing.")
  }

  finite_positive <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    is.finite(x) & x > 0
  }
  select_population_count <- function(tab) {
    corrected <- suppressWarnings(as.numeric(tab$correctedCount))
    raw_count <- suppressWarnings(as.numeric(tab$cellCount))
    countess <- suppressWarnings(as.numeric(tab$Countess))
    value <- ifelse(
      finite_positive(corrected), corrected,
      ifelse(
        finite_positive(raw_count), raw_count,
        ifelse(finite_positive(countess), countess, NA_real_)
      )
    )
    source <- ifelse(
      finite_positive(corrected), "correctedCount",
      ifelse(
        finite_positive(raw_count), "cellCount",
        ifelse(finite_positive(countess), "Countess", "missing")
      )
    )
    data.frame(value = value, source = source, stringsAsFactors = FALSE)
  }

  seed_rows <- population_snapshot[
    population_snapshot$event == "seeding" &
      population_snapshot$id %in% passage_ids,
    , drop = FALSE
  ]
  seed_selected <- select_population_count(seed_rows)
  seed_audit <- data.frame(
    passage_id = as.character(seed_rows$id),
    initial_population_size = seed_selected$value,
    initial_population_source = seed_selected$source,
    stringsAsFactors = FALSE
  )

  harvest_rows <- population_snapshot[
    population_snapshot$event == "harvest" &
      population_snapshot$passaged_from_id1 %in% passage_ids,
    , drop = FALSE
  ]
  harvest_selected <- select_population_count(harvest_rows)
  harvest_rows$selected_population_size <- harvest_selected$value
  harvest_rows$selected_population_source <- harvest_selected$source
  harvest_split <- split(
    harvest_rows,
    as.character(harvest_rows$passaged_from_id1),
    drop = TRUE
  )
  harvest_audit <- do.call(
    rbind,
    lapply(names(harvest_split), function(passage_id) {
      tab <- harvest_split[[passage_id]]
      data.frame(
        passage_id = passage_id,
        n_terminal_population_records = nrow(tab),
        terminal_population_size_median = stats::median(
          tab$selected_population_size, na.rm = TRUE
        ),
        terminal_population_size_min = min(
          tab$selected_population_size, na.rm = TRUE
        ),
        terminal_population_size_max = max(
          tab$selected_population_size, na.rm = TRUE
        ),
        terminal_population_sources = paste(
          sort(unique(tab$selected_population_source)),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(harvest_audit) <- NULL

  population_audit <- merge(
    passage[, c(
      "cohort", "segment_id", "passage_id",
      "cumulative_start_day", "cumulative_end_day"
    )],
    seed_audit,
    by = "passage_id", all.x = TRUE, sort = FALSE
  )
  population_audit <- merge(
    population_audit,
    harvest_audit,
    by = "passage_id", all.x = TRUE, sort = FALSE
  )
  population_audit <- population_audit[
    match(passage_ids, population_audit$passage_id),
    , drop = FALSE
  ]
  population_audit$initial_population_measured <- finite_positive(
    population_audit$initial_population_size
  )
  population_audit$terminal_population_measured <- finite_positive(
    population_audit$terminal_population_size_median
  )
  population_audit$population_measurement_complete <-
    population_audit$initial_population_measured &
    population_audit$terminal_population_measured
  population_audit$source_snapshot <- normalizePath(
    population_source, mustWork = TRUE
  )
  if (nrow(population_audit) != length(passage_ids) ||
      any(!population_audit$population_measurement_complete)) {
    stop(
      "Population-size audit did not confirm both initial and terminal counts ",
      "for all 131 Figure 1 in-vitro passage records."
    )
  }
  write_intermediate_tsv(
    population_audit,
    file.path(destination_root, "invitro_population_measurement_events.tsv")
  )

  required_flow_columns <- c(
    "segment_id", "cohort", "passage_id", "sample_name",
    "grid_index", "ploidy", "observed_density", "observed_log_density"
  )
  missing_flow_columns <- setdiff(required_flow_columns, names(observed_flow))
  if (length(missing_flow_columns)) {
    stop(
      "Observed flow table is missing columns: ",
      paste(missing_flow_columns, collapse = ", ")
    )
  }
  passage_match <- match(observed_flow$passage_id, passage$passage_id)
  if (anyNA(passage_match)) {
    stop("Observed flow samples did not all match Figure 1 passage records.")
  }
  observed_flow$cumulative_end_day <- passage$cumulative_end_day[passage_match]
  observed_flow$condition <- ifelse(
    observed_flow$oxygen_pct == 20.5,
    "control", "oxygen-deprived"
  )
  flow_split <- split(
    seq_len(nrow(observed_flow)), observed_flow$sample_name, drop = TRUE
  )
  observed_flow$normalized_density_mass <- NA_real_
  for (index in flow_split) {
    density <- suppressWarnings(as.numeric(observed_flow$observed_density[index]))
    if (any(!is.finite(density)) || sum(density) <= 0) {
      stop("Observed flow sample contains invalid density values.")
    }
    observed_flow$normalized_density_mass[index] <- density / sum(density)
  }
  if (length(flow_split) != 20L ||
      any(vapply(flow_split, length, integer(1L)) != 200L)) {
    stop("Expected 20 observed flow samples on a 200-point ploidy grid.")
  }

  flow_peak_rows <- do.call(
    rbind,
    lapply(flow_split, function(index) {
      tab <- observed_flow[index, , drop = FALSE]
      tab[which.max(tab$observed_density), c(
        "cohort", "condition", "segment_id", "passage_id",
        "sample_name", "cumulative_end_day", "ploidy", "observed_density"
      ), drop = FALSE]
    })
  )
  rownames(flow_peak_rows) <- NULL
  names(flow_peak_rows)[names(flow_peak_rows) == "ploidy"] <- "peak_ploidy"
  names(flow_peak_rows)[names(flow_peak_rows) == "observed_density"] <-
    "peak_observed_density"

  write_intermediate_tsv(
    observed_flow,
    file.path(destination_root, "invitro_flow_observed_density.tsv")
  )
  write_intermediate_tsv(
    flow_peak_rows,
    file.path(destination_root, "invitro_flow_peak_ploidy.tsv")
  )
  invisible(contract)
}

if (sys.nframe() == 0L) data_Figure1()
