#!/usr/bin/env Rscript

# Download the read-only CLONEID Passaging fields needed to construct the
# SUM-159 in-vitro dead-cell likelihood inputs. The resulting snapshots are
# deliberately self-contained so downstream work does not require DB access.

options(stringsAsFactors = FALSE, digits = 15)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve this script's path.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
output_dir <- dirname(script_path)
source_path <- file.path(output_dir, "sum159_all_db_cloneid_dead_cell_frequency.tsv")

if (!file.exists(source_path)) {
  stop("Missing dead-cell source TSV: ", source_path)
}

suppressPackageStartupMessages(library(cloneid))
suppressPackageStartupMessages(library(DBI))

dead <- read.delim(
  source_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = character()
)
if (!"db_source_id" %in% names(dead) || anyNA(dead$db_source_id) || any(!nzchar(dead$db_source_id))) {
  stop("Source TSV must contain non-empty db_source_id values.")
}
if (anyDuplicated(dead$db_source_id)) {
  stop("Source TSV contains duplicated db_source_id values.")
}

con <- cloneid::connect2DB()
on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

server <- DBI::dbGetQuery(
  con,
  paste(
    "SELECT DATABASE() AS database_name,",
    "NOW() AS server_time,",
    "@@session.time_zone AS session_time_zone,",
    "@@global.time_zone AS global_time_zone"
  )
)

required_passaging_columns <- c(
  "id", "cellLine", "event", "passaged_from_id1", "passaged_from_id2",
  "passage", "date"
)
optional_passaging_columns <- c(
  "cellCount", "correctedCount", "backup_cellCount", "Countess",
  "growthType", "comment", "media", "flask", "dishSurfaceArea_cm2"
)
available_passaging_columns <- DBI::dbListFields(con, "Passaging")
missing_required_columns <- setdiff(required_passaging_columns, available_passaging_columns)
if (length(missing_required_columns) > 0L) {
  stop(
    "Live Passaging table is missing required columns: ",
    paste(missing_required_columns, collapse = ", ")
  )
}
passaging_columns <- c(
  required_passaging_columns,
  intersect(optional_passaging_columns, available_passaging_columns)
)
quoted_columns <- paste(DBI::dbQuoteIdentifier(con, passaging_columns), collapse = ", ")
sum159_query <- paste0(
  "SELECT ", quoted_columns,
  " FROM Passaging WHERE cellLine = ", DBI::dbQuoteString(con, "SUM-159"),
  " ORDER BY date, id"
)
sum159 <- DBI::dbGetQuery(con, sum159_query)

quoted_ids <- paste(DBI::dbQuoteString(con, dead$db_source_id), collapse = ",")
source_query <- paste0(
  "SELECT ",
  "p.id AS db_source_id, p.cellLine AS db_cell_line, p.event AS db_event, ",
  "p.passage AS db_biological_passage, p.date AS db_event_datetime, ",
  "p.passaged_from_id1, p1.event AS parent1_event, p1.date AS parent1_datetime, ",
  "p.passaged_from_id2, p2.event AS parent2_event, p2.date AS parent2_datetime, ",
  "TIMESTAMPDIFF(SECOND, p1.date, p.date) / 86400.0 AS elapsed_days_from_parent1, ",
  "TIMESTAMPDIFF(SECOND, p2.date, p.date) / 86400.0 AS elapsed_days_from_parent2 ",
  "FROM Passaging p ",
  "LEFT JOIN Passaging p1 ON p1.id = p.passaged_from_id1 ",
  "LEFT JOIN Passaging p2 ON p2.id = p.passaged_from_id2 ",
  "WHERE p.id IN (", quoted_ids, ") ",
  "ORDER BY p.date, p.id"
)
source_rows <- DBI::dbGetQuery(con, source_query)

if (nrow(sum159) == 0L) {
  stop("The SUM-159 Passaging query returned no rows.")
}
if (nrow(source_rows) != nrow(dead) || !setequal(source_rows$db_source_id, dead$db_source_id)) {
  missing_ids <- setdiff(dead$db_source_id, source_rows$db_source_id)
  extra_ids <- setdiff(source_rows$db_source_id, dead$db_source_id)
  stop(
    "The source-specific DB snapshot did not match the source TSV. Missing: ",
    paste(missing_ids, collapse = ", "), "; extra: ", paste(extra_ids, collapse = ", ")
  )
}
if (!all(dead$db_source_id %in% sum159$id)) {
  stop(
    "Some dead-cell db_source_id values are absent from the complete SUM-159 snapshot: ",
    paste(setdiff(dead$db_source_id, sum159$id), collapse = ", ")
  )
}

source_rows <- source_rows[match(dead$db_source_id, source_rows$db_source_id), , drop = FALSE]
if (!identical(source_rows$db_source_id, dead$db_source_id)) {
  stop("Failed to restore source TSV row order after the DB query.")
}
if (!identical(as.character(source_rows$db_event), as.character(dead$event))) {
  bad <- which(as.character(source_rows$db_event) != as.character(dead$event))
  stop("DB/source event mismatch for: ", paste(dead$db_source_id[bad], collapse = ", "))
}

stamp <- format(as.POSIXct(server$server_time[[1]], tz = "UTC"), "%Y%m%d")
sum159_path <- file.path(output_dir, paste0("cloneid_passaging_sum159_snapshot_", stamp, ".tsv"))
source_rows_path <- file.path(
  output_dir,
  paste0("cloneid_passaging_dead_cell_sources_and_parents_", stamp, ".tsv")
)
metadata_path <- file.path(output_dir, paste0("cloneid_passaging_snapshot_metadata_", stamp, ".tsv"))

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

write_tsv(sum159, sum159_path)
write_tsv(source_rows, source_rows_path)

metadata <- data.frame(
  key = c(
    "database_name", "database_server_time", "database_session_time_zone",
    "database_global_time_zone", "source_dead_cell_tsv", "source_dead_cell_tsv_md5",
    "sum159_snapshot_file", "sum159_snapshot_rows", "sum159_snapshot_min_datetime",
    "sum159_snapshot_max_datetime", "source_parent_snapshot_file",
    "source_parent_snapshot_rows", "source_ids_expected", "source_ids_matched",
    "sum159_query", "source_parent_query"
  ),
  value = c(
    as.character(server$database_name[[1]]),
    as.character(server$server_time[[1]]),
    as.character(server$session_time_zone[[1]]),
    as.character(server$global_time_zone[[1]]),
    basename(source_path),
    unname(tools::md5sum(source_path)),
    basename(sum159_path),
    nrow(sum159),
    as.character(min(sum159$date)),
    as.character(max(sum159$date)),
    basename(source_rows_path),
    nrow(source_rows),
    nrow(dead),
    sum(dead$db_source_id %in% source_rows$db_source_id),
    sum159_query,
    source_query
  ),
  stringsAsFactors = FALSE
)
write_tsv(metadata, metadata_path)

cat("SUM159_SNAPSHOT=", sum159_path, "\n", sep = "")
cat("SOURCE_PARENT_SNAPSHOT=", source_rows_path, "\n", sep = "")
cat("METADATA=", metadata_path, "\n", sep = "")
cat("SUM159_ROWS=", nrow(sum159), "\n", sep = "")
cat("SOURCE_ROWS=", nrow(source_rows), "\n", sep = "")
