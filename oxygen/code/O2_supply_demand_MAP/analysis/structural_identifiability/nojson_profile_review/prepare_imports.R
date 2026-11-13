# Numerical artifact import only. No simulator or optimizer is sourced.
args <- commandArgs(TRUE); stopifnot(length(args) == 2L)
project <- normalizePath(args[1]); study <- normalizePath(args[2])
tab <- function(x) read.delim(x, stringsAsFactors = FALSE, check.names = FALSE)
put <- function(x, path) {
  x[] <- lapply(x, function(v) if (is.double(v)) {z <- sprintf("%.17g", v); z[is.na(v)] <- NA_character_; z} else v)
  write.table(x, path, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")
}
clean_characters <- function(x) {
  if (is.character(x)) {
    z <- x[!is.na(x)]
    stopifnot(!any(grepl("^[[:space:]]*(\\{[[:space:]]*\"|\\[[[:space:]]*[\\{\"0-9-])", z)),
              !any(grepl("[\u3400-\u4dbf\u4e00-\u9fff]", z)))
  } else if (is.list(x)) invisible(lapply(x, clean_characters))
  invisible(TRUE)
}
review <- "oxygen/results/structural_identifiability/20261110_culture_v2_optimization_review"
baseline <- "oxygen/results/structural_identifiability/20261111_culture_v2_baseline_recovery_nojson"
q <- tab(file.path(project, review, "tables/review_point_quality.tsv"))
order <- c("pm_00", "bb_00", "pm_05", "bb_07", "pm_06", "bb_08", "pm_09")
q <- q[match(order, q$point_id), ]; stopifnot(!anyNA(q$point_id))
registry <- data.frame(point = q$point_id, target = q$target, value = q$value,
                       priority = seq_len(nrow(q)), historical_task = q$best_task,
                       historical_objective = q$objective)
sources <- data.frame(import_id = paste0("historical_", q$point_id), point = q$point_id,
                      origin = q$best_task, family = "", kind = "historical_point",
                      source = file.path(review, "fits", q$best_task, "accepted_fit.rds"))
ids <- c("base_reference", "base_A1", "base_A2", "base_B2", "base_B4")
fits <- c("consistency_1", "A1_local", "consistency_2", "consistency_3", "consistency_4")
sources <- rbind(sources, data.frame(import_id = ids, point = "base", origin = c("review_035", "A1", "A2", "B2", "B4"),
  family = c("", "A", "A", "B", "B"), kind = c("historical_reference", rep("imported_independent_baseline", 4)),
  source = file.path(baseline, "fits", fits, "accepted_fit.rds")))
dir.create(file.path(study, "imports"), showWarnings = FALSE)
for (i in seq_len(nrow(sources))) {
  r <- sources[i, ]; fit <- readRDS(file.path(project, r$source))
  if (r$point != "base") {
    j <- match(r$point, registry$point)
    stopifnot(identical(registry$value[j], fit$best_params[[registry$target[j]]]),
              abs(registry$historical_objective[j] - fit$best_components$objective) <= 1e-12)
  }
  imported <- list(best_params = fit$best_params, best_components = fit$best_components)
  clean_characters(imported)
  sources$objective[i] <- fit$best_components$objective
  sources$relative_path[i] <- file.path("imports", paste0(r$import_id, ".rds"))
  saveRDS(imported, file.path(study, sources$relative_path[i]))
}
put(registry, file.path(study, "cases.tsv")); put(sources, file.path(study, "imports.tsv"))
cat("Imported twelve clean numerical RDS artifacts; seven fixed targets exactly match original numerical values.\n")
