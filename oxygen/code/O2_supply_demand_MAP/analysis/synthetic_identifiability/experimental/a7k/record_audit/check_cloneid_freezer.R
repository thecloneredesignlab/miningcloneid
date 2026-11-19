# Read only the donor, revival, and harvest records in scope.
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1L]]
cohorts <- c("2N_O1", "2N_O2", "4N_O1", "4N_O2")
donors <- paste0("SUM-159_NLS_", cohorts, "_A7_seedT2")
seeds <- paste0("SUM-159_NLS_", cohorts, "_A7K_seed")
a8 <- paste0("SUM-159_NLS_", cohorts, "_A8_seed")
harvests <- paste0("SUM159_NLS_", cohorts, "_A7K_harvest")
export <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE,
              quote = FALSE, na = "")
}
run <- function() {
  suppressPackageStartupMessages(suppressWarnings(requireNamespace("cloneid")))
  detail <- "connection_error"
  con <- tryCatch(suppressPackageStartupMessages(cloneid::connect2DB()), error=function(e) {
    if (grepl("access denied", conditionMessage(e), ignore.case=TRUE)) {
      detail <<- "access_denied"
    }
    NULL
  })
  if (is.null(con)) {
    export(data.frame(status="configured_connection_unavailable", detail=detail),
           "cloneid_query_status.tsv")
    return(invisible(NULL))
  }
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tables <- DBI::dbGetQuery(con, paste0(
    "SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA=DATABASE() ",
    "AND LOWER(TABLE_NAME) REGEXP 'freez|cry|thaw|vial|storage' ORDER BY TABLE_NAME"))
  export(tables, "cloneid_storage_tables.tsv")
  in_sql <- function(ids) paste(DBI::dbQuoteString(con, ids), collapse=",")
  ids <- c(donors, seeds, a8, harvests)
  if ("Minus80Freezer" %in% tables$TABLE_NAME) {
    freeze <- DBI::dbGetQuery(con, paste0(
      "SELECT id, cellCount, Comment FROM Minus80Freezer WHERE id IN (", in_sql(ids),
      ") ORDER BY id"))
    export(freeze,"cloneid_freezer_exact_ids.tsv")
    cat("Exact freezer ledger rows", nrow(freeze), "\n")
  }
  fields <- DBI::dbListFields(con, "Passaging")
  wanted <- c("id","event","date","passaged_from_id1","passaged_from_id2",
              "media","cellCount","correctedCount","comment")
  extras <- grep("freez|cry|thaw|vial|storage", fields, value=TRUE, ignore.case=TRUE)
  wanted <- unique(c(intersect(wanted,fields), extras))
  col_sql <- paste(DBI::dbQuoteIdentifier(con,wanted), collapse=",")
  invariant_ids <- c(seeds, a8, harvests)
  where <- paste0(" WHERE id IN (",in_sql(ids),") OR passaged_from_id1 IN (",
                  in_sql(donors),") OR passaged_from_id2 IN (",in_sql(donors),")")
  events <- DBI::dbGetQuery(con,paste0("SELECT ",col_sql," FROM Passaging",where," ORDER BY date,id"))
  export(events,"cloneid_live_lineage_events.tsv")
  cat("Exact lineage events", nrow(events), "; fields about storage in Passaging",length(extras),"\n")
  cat("Related tables",paste(tables$TABLE_NAME,collapse=", "), "\n")
  export(data.frame(status="read_only_queries_completed",exact_anchor_rows=sum(events$id %in% ids),
                    children_other_than_registered=sum(!events$id %in% c(donors,invariant_ids))),
         "cloneid_query_status.tsv")
}
tryCatch(run(),error=function(e) {
  export(data.frame(status="query_unavailable",error_category=class(e)[1]), "cloneid_query_status.tsv")
  stop("Read-only CloneID query unavailable; see query status.", call.=FALSE)
})
