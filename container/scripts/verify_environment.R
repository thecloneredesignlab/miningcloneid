#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: verify_environment.R PACKAGES_LOCK_TSV [hpc-observed|target]", call. = FALSE)
}

lock_path <- args[[1L]]
profile <- if (length(args) >= 2L) args[[2L]] else "target"
if (!profile %in% c("hpc-observed", "target")) {
  stop("Profile must be hpc-observed or target", call. = FALSE)
}

lock <- utils::read.delim(lock_path, check.names = FALSE, stringsAsFactors = FALSE)
expected_r <- "4.4.2"
r_ok <- identical(paste(R.version$major, R.version$minor, sep = "."), expected_r)

rows <- lapply(seq_len(nrow(lock)), function(i) {
  package <- lock$package[[i]]
  expected <- if (profile == "hpc-observed") {
    lock$observed_hpc_version[[i]]
  } else {
    lock$target_version[[i]]
  }
  if (is.na(expected) || !nzchar(expected)) return(NULL)
  installed <- requireNamespace(package, quietly = TRUE)
  actual <- if (installed) {
    as.character(utils::packageDescription(package, fields = "Version"))
  } else {
    NA_character_
  }
  data.frame(
    package = package,
    expected = expected,
    actual = actual,
    installed = installed,
    version_match = installed && identical(actual, expected),
    stringsAsFactors = FALSE
  )
})
rows <- Filter(Negate(is.null), rows)
result <- do.call(rbind, rows)
utils::write.table(
  result,
  file = stdout(),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)

if (!r_ok || any(!result$version_match)) {
  quit(status = 1L, save = "no")
}
