#!/usr/bin/env Rscript

# Canonical fitted-result discovery, parameter input, and model-loading helpers.

o2ipa_find_extra <- function(run_dir, file) {
    path <- file.path(run_dir, "extra_results", file)
    if (file.exists(path))
        path
    else NA_character_
}


o2ipa_read_extra_summary <- function(run_dir) {
    path <- o2ipa_find_extra(run_dir, "seed_summary.tsv")
    if (is.na(path))
        return(NULL)
    tab <- o2ipa_read_tsv(path)
    if (!"seed" %in% names(tab))
        return(NULL)
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    tab
}


o2ipa_read_param_matrix <- function(run_dir) {
    parent <- dirname(normalizePath(run_dir, mustWork = FALSE))
    base <- basename(normalizePath(run_dir, mustWork = FALSE))
    candidates <- c(file.path(run_dir, "param_matrix_with_seed.tsv"), file.path(run_dir, "parameter_matrix_with_seed.tsv"),
        file.path(parent, paste0(base, "_param_matrix_with_seed.tsv")), file.path(parent, paste0(base, "_param_matrix.tsv")))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates))
        return(NULL)
    tab <- o2ipa_read_tsv(candidates[[1]])
    if (!"seed" %in% names(tab))
        tab$seed <- seq_len(nrow(tab))
    tab$seed_id <- o2ipa_norm_seed(tab$seed)
    attr(tab, "source_file") <- candidates[[1]]
    tab
}


o2ipa_discover_seeds <- function(run_dir) {
    dirs <- list.dirs(run_dir, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[grepl("^seed[^/]*[0-9]+$", basename(dirs))]
    seed_ids <- o2ipa_norm_seed(basename(dirs))
    seed_dir_map <- setNames(normalizePath(dirs, mustWork = FALSE), seed_ids)
    summary_tab <- o2ipa_read_extra_summary(run_dir)
    matrix_tab <- o2ipa_read_param_matrix(run_dir)
    all_ids <- unique(c(seed_ids, if (!is.null(summary_tab)) summary_tab$seed_id, if (!is.null(matrix_tab)) matrix_tab$seed_id))
    all_ids <- all_ids[o2ipa_order_seeds(all_ids)]
    data.frame(seed_id = all_ids, seed_dir = unname(seed_dir_map[all_ids]), stringsAsFactors = FALSE)
}


o2ipa_metric_map <- function(path) {
    if (is.null(path) || is.na(path) || !file.exists(path))
        return(list())
    tab <- o2ipa_read_tsv(path)
    if (!all(c("metric", "value") %in% names(tab)))
        return(list())
    vals <- as.list(tab$value)
    names(vals) <- tab$metric
    vals
}


o2ipa_num_from_map <- function(map, key) {
    if (is.null(map[[key]]))
        return(NA_real_)
    suppressWarnings(as.numeric(map[[key]]))
}


o2ipa_chr_from_map <- function(map, key) {
    if (is.null(map[[key]]))
        return(NA_character_)
    as.character(map[[key]])
}


o2ipa_read_best_params <- function(path) {
    if (!file.exists(path))
        return(list(values = setNames(numeric(0), character(0)), aliases = data.frame()))
    tab <- o2ipa_read_tsv(path)
    if (!all(c("parameter", "value") %in% names(tab))) {
        stop("best_params.tsv must contain parameter and value columns: ", path)
    }
    vals <- suppressWarnings(as.numeric(tab$value))
    names(vals) <- trimws(as.character(tab$parameter))
    list(values = vals, aliases = data.frame())
}


o2ipa_extract_param <- function(seed_id, parameter, best_vals, summary_row, matrix_row) {
    aliases <- o2ipa_param_aliases()[[parameter]]
    sources <- list()
    if (length(best_vals)) {
        sources$best_params <- best_vals
    }
    if (!is.null(summary_row) && nrow(summary_row) == 1L) {
        vals <- numeric(0)
        for (a in aliases) {
            col <- paste0("value__", a)
            if (col %in% names(summary_row))
                vals[a] <- suppressWarnings(as.numeric(summary_row[[col]][[1]]))
        }
        sources$seed_summary_value_cols <- vals
    }
    if (!is.null(matrix_row) && nrow(matrix_row) == 1L) {
        vals <- numeric(0)
        for (a in aliases) {
            for (col in c(paste0("final_", a), a)) {
                if (col %in% names(matrix_row))
                  vals[a] <- suppressWarnings(as.numeric(matrix_row[[col]][[1]]))
            }
        }
        sources$param_matrix <- vals
    }
    for (src in names(sources)) {
        vals <- sources[[src]]
        for (a in aliases) {
            if (a %in% names(vals)) {
                v <- suppressWarnings(as.numeric(vals[[a]]))
                if (is.finite(v)) {
                  return(list(value = v, source = src, alias = a))
                }
            }
        }
    }
    list(value = NA_real_, source = NA_character_, alias = NA_character_)
}


o2ipa_choose_objective <- function(summary_vals, objective_source = "auto") {
    raw <- NA_real_
    b <- o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw")
    p <- o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw")
    if (is.finite(b) && is.finite(p))
        raw <- b + p
    data <- o2ipa_num_from_map(summary_vals, "objective_data")
    total <- o2ipa_num_from_map(summary_vals, "objective_total")
    if (!is.finite(total))
        total <- o2ipa_num_from_map(summary_vals, "objective")
    source <- match.arg(objective_source, c("auto", "raw_likelihood", "data", "total"))
    if (identical(source, "raw_likelihood"))
        return(list(value = raw, source = "raw_likelihood"))
    if (identical(source, "data"))
        return(list(value = data, source = "objective_data"))
    if (identical(source, "total"))
        return(list(value = total, source = if (is.finite(o2ipa_num_from_map(summary_vals, "objective_total"))) "objective_total" else "objective"))
    if (is.finite(raw))
        return(list(value = raw, source = "raw_likelihood"))
    if (is.finite(data))
        return(list(value = data, source = "objective_data"))
    list(value = total, source = if (is.finite(total)) "objective_total_or_objective" else NA_character_)
}


o2ipa_source_model <- function(script_dir = o2ipa_script_dir()) {
    workflow_root <- o2ipa_workflow_root(script_dir)
    model_file <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
    if (!file.exists(model_file))
        stop("Model file not found: ", model_file)
    shared_file <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
    common_file <- file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R")
    if (file.exists(shared_file))
        source(shared_file, local = globalenv())
    if (file.exists(common_file))
        source(common_file, local = globalenv())
    env <- new.env(parent = globalenv())
    old_dir <- Sys.getenv("MININGCLONEID_OXYGEN_CODE_DIR", unset = NA_character_)
    Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_file))
    on.exit({
        if (is.na(old_dir)) Sys.unsetenv("MININGCLONEID_OXYGEN_CODE_DIR") else Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = old_dir)
    }, add = TRUE)
    source(model_file, local = env)
    env
}


o2ipa_call_model <- function(model_env, fn, ...) {
    if (!exists(fn, envir = model_env, inherits = TRUE))
        stop("Model helper missing: ", fn)
    get(fn, envir = model_env, inherits = TRUE)(...)
}
