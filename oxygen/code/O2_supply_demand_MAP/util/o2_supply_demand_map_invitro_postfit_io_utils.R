# Shared post-fit in-vitro I/O, payload extraction, lineage metadata, and
# schema helpers used across simulation domains.
# This module does not plot or execute simulations.

ivt_sim_first_non_null <- function(...) {
  values <- list(...)
  for (value in values) {
    if (!is.null(value) && length(value)) return(value)
  }
  NULL
}

ivt_sim_write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  invisible(path)
}

ivt_sim_extract_payload <- function(fit_result, fit_dir, best_params_path = NULL) {
  if (!is.list(fit_result)) stop("fit_result.rds must contain a list.")
  best_components <- fit_result$best_components
  comp <- best_components
  if (is.list(comp) && !is.null(comp$invitro)) comp <- comp$invitro
  if (!is.list(comp)) stop("fit_result.rds does not contain in-vitro best_components.")

  cfg <- ivt_sim_first_non_null(
    fit_result$cfg,
    fit_result$ctx$invitro_cfg,
    fit_result$ctx$invitro$cfg,
    comp$cfg
  )
  if (!is.list(cfg)) stop("fit_result.rds does not contain an in-vitro cfg.")
  cfg <- normalize_sim_cfg_common(cfg, context = "viz")

  run_params <- ivt_sim_first_non_null(
    fit_result$best_params,
    fit_result$invitro_run_params,
    best_components$invitro_run_params,
    comp$invitro_run_params,
    comp$run_params,
    fit_result$ctx$invitro_run_params,
    fit_result$ctx$invitro$run_params
  )

  if (is.null(best_params_path)) {
    candidate <- file.path(fit_dir, "best_params.tsv")
    if (file.exists(candidate)) best_params_path <- candidate
  }
  if (!is.null(best_params_path) && file.exists(best_params_path)) {
    best_df <- utils::read.delim(
      best_params_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (all(c("parameter", "value") %in% names(best_df))) {
      if (!is.list(run_params)) run_params <- list()
      for (i in seq_len(nrow(best_df))) {
        name <- as.character(best_df$parameter[[i]])
        value <- suppressWarnings(as.numeric(best_df$value[[i]]))
        if (nzchar(name) && is.finite(value)) run_params[[name]] <- value
      }
    }
  }
  if (!is.list(run_params)) stop("fit_result.rds does not contain in-vitro best parameters.")
  run_params <- normalize_run_params_common(run_params, cfg = cfg)

  list(
    cfg = cfg,
    run_params = run_params,
    components = comp,
    best_params_path = best_params_path
  )
}

ivt_sim_extract_optimizer_population <- function(fit_result) {
  candidates <- list(
    fit_result$deoptim$member$bestmemit,
    fit_result$deoptim$member$pop,
    fit_result$optimizer_trace
  )
  for (candidate in candidates) {
    if (is.null(candidate) || !(is.matrix(candidate) || is.data.frame(candidate))) next
    population <- as.data.frame(candidate, stringsAsFactors = FALSE)
    population[] <- lapply(population, function(x) suppressWarnings(as.numeric(x)))
    keep <- vapply(
      population,
      function(x) sum(is.finite(x)) >= 3L && stats::var(x, na.rm = TRUE) > 0,
      logical(1)
    )
    population <- population[, keep, drop = FALSE]
    population <- population[stats::complete.cases(population), , drop = FALSE]
    if (nrow(population) >= 3L && ncol(population) >= 2L) {
      return(data.frame(
        optimizer_sample_id = seq_len(nrow(population)),
        population,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ))
    }
  }
  data.frame(optimizer_sample_id = integer(), stringsAsFactors = FALSE)
}

ivt_sim_lineage_label <- function(key) {
  key <- as.character(key)
  out <- rep(NA_character_, length(key))
  ok <- !is.na(key) & nzchar(key)
  if (!any(ok)) return(out)
  parts <- strsplit(key[ok], "_", fixed = TRUE)
  control <- vapply(
    parts,
    function(x) length(x) > 0L && all(trimws(x) == "20.5"),
    logical(1)
  )
  out[ok] <- ifelse(control, "control", "deprived")
  out
}

ivt_sim_normalize_lineage_columns <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  n <- nrow(df)
  if (!"segment_id" %in% names(df)) df$segment_id <- as.character(seq_len(n))
  if (!"lineage_terminal_key" %in% names(df)) {
    df$lineage_terminal_key <- as.character(df$segment_id)
  }
  if (!"lineage_passage_index" %in% names(df)) {
    df$lineage_passage_index <- if ("passage_index" %in% names(df)) {
      suppressWarnings(as.numeric(df$passage_index))
    } else {
      seq_len(n)
    }
  }
  if (!"lineage_label" %in% names(df)) {
    df$lineage_label <- ivt_sim_lineage_label(df$lineage_terminal_key)
  }
  missing_label <- is.na(df$lineage_label) | !nzchar(as.character(df$lineage_label))
  if (any(missing_label)) {
    fallback <- if ("cohort" %in% names(df)) as.character(df$cohort) else rep("lineage", n)
    df$lineage_label[missing_label] <- fallback[missing_label]
  }
  if (!"parent_segment_id" %in% names(df)) df$parent_segment_id <- NA_character_
  if (!"selected_day" %in% names(df)) df$selected_day <- NA_real_
  df
}

ivt_sim_schema <- function(tables) {
  dplyr::bind_rows(lapply(names(tables), function(table_name) {
    table <- tables[[table_name]]
    data.frame(
      table = paste0(table_name, ".tsv"),
      column = names(table),
      class = vapply(table, function(x) paste(class(x), collapse = ","), character(1)),
      rows = nrow(table),
      stringsAsFactors = FALSE
    )
  }))
}
