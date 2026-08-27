#!/usr/bin/env Rscript

# Shared fixed-O2 fitting-parameter and model-loading helpers used by
# simulation/fix_o2_simulation.R and analyses that source it.
suppressPackageStartupMessages({
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
})

o2ipa_null_coalesce <- function(x, y) {
    if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
        y
    else x
}

o2ipa_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]
        if (!grepl("^--", arg)) {
            i <- i + 1L
            next
        }
        kv <- sub("^--", "", arg)
        eq <- regexpr("=", kv, fixed = TRUE)
        if (eq > 0L) {
            key <- substr(kv, 1L, eq - 1L)
            val <- substr(kv, eq + 1L, nchar(kv))
            out[[key]] <- val
            i <- i + 1L
        } else {
            key <- kv
            if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
                out[[key]] <- args[[i + 1L]]
                i <- i + 2L
            } else {
                out[[key]] <- TRUE
                i <- i + 1L
            }
        }
    }
    out
}

o2ipa_as_chr <- function(x, default = "") {
    val <- o2ipa_null_coalesce(x, default)
    val <- as.character(val[[1]])
    if (!nzchar(val))
        default
    else val
}

o2ipa_as_num <- function(x, default = NA_real_) {
    val <- suppressWarnings(as.numeric(o2ipa_null_coalesce(x, default)[[1]]))
    if (is.finite(val))
        val
    else default
}

o2ipa_as_int <- function(x, default = NA_integer_) {
    val <- suppressWarnings(as.integer(o2ipa_null_coalesce(x, default)[[1]]))
    if (!is.na(val))
        val
    else default
}

o2ipa_as_bool <- function(x, default = FALSE) {
    if (is.null(x) || !length(x) || is.na(x[[1]]))
        return(isTRUE(default))
    if (is.logical(x[[1]]))
        return(isTRUE(x[[1]]))
    tolower(trimws(as.character(x[[1]]))) %in% c("1", "true", "t", "yes", "y", "on")
}

o2ipa_split_csv <- function(x, default = character()) {
    txt <- trimws(o2ipa_as_chr(x, paste(default, collapse = ",")))
    if (!nzchar(txt))
        return(default)
    vals <- trimws(strsplit(txt, ",", fixed = TRUE)[[1]])
    vals[nzchar(vals)]
}

o2ipa_script_dir <- function() {
    frame_files <- Filter(
        nzchar,
        vapply(sys.frames(), function(env) {
            ofile <- env$ofile
            if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
        }, character(1))
    )
    if (length(frame_files)) {
        return(dirname(frame_files[[length(frame_files)]]))
    }
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)))
    }
    normalizePath(getwd(), mustWork = FALSE)
}

o2ipa_find_workflow_root <- function(path = o2ipa_script_dir()) {
    cur <- normalizePath(path, mustWork = FALSE)
    if (file.exists(cur) && !dir.exists(cur))
        cur <- dirname(cur)
    for (i in seq_len(8L)) {
        if (file.exists(file.path(cur, "util", "o2_supply_demand_map_shared.R")) && file.exists(file.path(cur, "model", "model_O2_supply_demand_MAP.R"))) {
            return(normalizePath(cur, mustWork = FALSE))
        }
        repo_candidate <- file.path(cur, "oxygen", "code", "O2_supply_demand_MAP")
        if (file.exists(file.path(repo_candidate, "util", "o2_supply_demand_map_shared.R")) && file.exists(file.path(repo_candidate, "model", "model_O2_supply_demand_MAP.R"))) {
            return(normalizePath(repo_candidate, mustWork = FALSE))
        }
        parent <- dirname(cur)
        if (identical(parent, cur))
            break
        cur <- parent
    }
    normalizePath(file.path(path, "..", ".."), mustWork = FALSE)
}

o2ipa_workflow_root <- function(script_dir = o2ipa_script_dir()) {
    o2ipa_find_workflow_root(script_dir)
}

o2ipa_repo_root <- function(script_dir = o2ipa_script_dir()) {
    normalizePath(file.path(o2ipa_workflow_root(script_dir), "..", "..", ".."), mustWork = FALSE)
}

o2ipa_read_tsv <- function(path) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
}

o2ipa_seed_number <- function(seed_id) {
    suppressWarnings(as.integer(sub("^seed", "", as.character(seed_id))))
}

o2ipa_norm_seed <- function(x) {
    s <- as.character(x)
    s <- trimws(s)
    out <- ifelse(grepl("^seed[0-9]+$", s), s, ifelse(grepl("^[0-9]+$", s), paste0("seed", s), s))
    out
}

o2ipa_order_seeds <- function(seed_id) {
    n <- o2ipa_seed_number(seed_id)
    order(ifelse(is.na(n), Inf, n), seed_id)
}

o2ipa_target_params <- function() {
    c("O2_crit", "alpha_o2", "gamma_growth", "mu_hp", "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta", "buffer_n_exp",
        "n_O", "lam_max", "p_mis_base", "p_wgd", "gamma_mu", "o2_S0", "kappa_O", "eta_o2", "rho_2N")
}

o2ipa_param_aliases <- function() {
    list(O2_crit = c("O2_crit", "o2_crit", "O2crit", "o2crit"), alpha_o2 = c("alpha_o2", "alpha_O2"), gamma_growth = c("gamma_growth"),
        mu_hp = c("mu_hp"), p_misseg = c("p_misseg", "p_mis", "p_missegregation"), k_o_mis = c("k_o_mis", "ko_mis", "k_O_mis"),
        buffer_smax = c("buffer_smax", "buffer_s_max"), buffer_beta = c("buffer_beta"), buffer_n_exp = c("buffer_n_exp",
            "buffer_n"), n_O = c("n_O", "n_o", "nO"), lam_max = c("lam_max", "lambda_max"), p_mis_base = c("p_mis_base",
            "p_misseg_base"), p_wgd = c("p_wgd", "wgd_prob"), gamma_mu = c("gamma_mu"), o2_S0 = c("o2_S0", "O2_S0", "S0_o2"),
        kappa_O = c("kappa_O", "kappa_o"), eta_o2 = c("eta_o2", "eta_O2"), rho_2N = c("rho_2N", "rho2N"))
}

o2ipa_optimizer_transform <- function(parameter) {
    log10_params <- c("O2_crit", "alpha_o2", "mu_hp", "p_misseg", "k_o_mis", "buffer_beta", "buffer_n_exp", "lam_max", "p_mis_base",
        "p_wgd", "o2_S0", "kappa_O", "eta_o2", "rho_2N")
    if (parameter %in% log10_params)
        "log10"
    else "identity"
}

o2ipa_transform_parameter_value <- function(parameter, value, epsilon = 1e-12) {
    v <- as.numeric(value)
    if (!is.finite(v))
        return(NA_real_)
    tr <- o2ipa_optimizer_transform(parameter)
    if (identical(tr, "log10")) {
        if (v <= 0)
            return(NA_real_)
        return(log10(v))
    }
    if (identical(tr, "logit")) {
        vv <- min(max(v, epsilon), 1 - epsilon)
        return(stats::qlogis(vv))
    }
    v
}

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

o2ipa_extract_all_params <- function(seed_id, seed_dir, summary_tab, matrix_tab) {
    best_path <- if (!is.na(seed_dir) && nzchar(seed_dir))
        file.path(seed_dir, "best_params.tsv")
    else NA_character_
    best_vals <- if (!is.na(best_path) && file.exists(best_path))
        o2ipa_read_best_params(best_path)$values
    else setNames(numeric(0), character(0))
    summary_row <- if (!is.null(summary_tab))
        summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE]
    else NULL
    if (!is.null(summary_row) && nrow(summary_row) > 1L)
        summary_row <- summary_row[1, , drop = FALSE]
    matrix_row <- if (!is.null(matrix_tab))
        matrix_tab[matrix_tab$seed_id == seed_id, , drop = FALSE]
    else NULL
    if (!is.null(matrix_row) && nrow(matrix_row) > 1L)
        matrix_row <- matrix_row[1, , drop = FALSE]
    params <- o2ipa_target_params()
    rows <- lapply(params, function(p) {
        hit <- o2ipa_extract_param(seed_id, p, best_vals, summary_row, matrix_row)
        data.frame(seed_id = seed_id, parameter = p, value = hit$value, parameter_source = hit$source, matched_alias = hit$alias,
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
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

o2ipa_collect_seed_inputs <- function(run_dir, objective_source = "auto") {
    manifest0 <- o2ipa_discover_seeds(run_dir)
    summary_tab <- o2ipa_read_extra_summary(run_dir)
    matrix_tab <- o2ipa_read_param_matrix(run_dir)
    boundary_path <- o2ipa_find_extra(run_dir, "parameter_boundary_long.tsv")
    boundary_long <- if (!is.na(boundary_path))
        o2ipa_read_tsv(boundary_path)
    else NULL
    if (!is.null(boundary_long) && "seed" %in% names(boundary_long)) {
        boundary_long$seed_id <- o2ipa_norm_seed(boundary_long$seed)
    }
    param_rows <- list()
    manifest_rows <- vector("list", nrow(manifest0))
    for (i in seq_len(nrow(manifest0))) {
        seed_id <- manifest0$seed_id[[i]]
        seed_dir <- manifest0$seed_dir[[i]]
        if (is.na(seed_dir))
            seed_dir <- ""
        fit_summary_path <- if (nzchar(seed_dir))
            file.path(seed_dir, "fit_summary.tsv")
        else NA_character_
        summary_vals <- if (!is.na(fit_summary_path) && file.exists(fit_summary_path)) {
            o2ipa_metric_map(fit_summary_path)
        }
        else {
            row <- if (!is.null(summary_tab))
                summary_tab[summary_tab$seed_id == seed_id, , drop = FALSE]
            else NULL
            if (!is.null(row) && nrow(row))
                as.list(row[1, , drop = TRUE])
            else list()
        }
        obj <- o2ipa_choose_objective(summary_vals, objective_source = objective_source)
        params_long <- o2ipa_extract_all_params(seed_id, seed_dir, summary_tab, matrix_tab)
        param_rows[[i]] <- params_long
        missing_params <- params_long$parameter[!is.finite(params_long$value)]
        bseed <- if (!is.null(boundary_long))
            boundary_long[boundary_long$seed_id == seed_id, , drop = FALSE]
        else NULL
        target_boundary <- if (!is.null(bseed) && nrow(bseed)) {
            pname <- if ("param_prototype" %in% names(bseed))
                bseed$param_prototype
            else if ("parameter" %in% names(bseed))
                bseed$parameter
            else bseed$param_name
            bseed[pname %in% o2ipa_target_params(), , drop = FALSE]
        }
        else {
            data.frame()
        }
        boundary_risk <- FALSE
        n_near <- 0L
        if (nrow(target_boundary)) {
            status_col <- if ("bound_status" %in% names(target_boundary))
                "bound_status"
            else if ("status" %in% names(target_boundary))
                "status"
            else NA_character_
            if (!is.na(status_col)) {
                boundary_risk <- any(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior",
                  ""))
                n_near <- sum(!is.na(target_boundary[[status_col]]) & !target_boundary[[status_col]] %in% c("interior", ""))
            }
        }
        fit_success <- length(missing_params) == 0L && is.finite(obj$value)
        convergence <- o2ipa_chr_from_map(summary_vals, "deoptim_stop_reason")
        if (is.na(convergence))
            convergence <- o2ipa_chr_from_map(summary_vals, "optimizer_local_convergence")
        failure_parts <- character(0)
        if (!is.finite(obj$value))
            failure_parts <- c(failure_parts, "missing_objective")
        if (length(missing_params))
            failure_parts <- c(failure_parts, paste0("missing_params:", paste(missing_params, collapse = ",")))
        manifest_rows[[i]] <- data.frame(seed_id = seed_id, seed_dir = if (nzchar(seed_dir))
            normalizePath(seed_dir, mustWork = FALSE)
        else NA_character_, fit_success = fit_success, convergence_status = convergence, objective = obj$value, objective_source = obj$source,
            objective_total = o2ipa_num_from_map(summary_vals, "objective_total"), objective_data = o2ipa_num_from_map(summary_vals,
                "objective_data"), objective_burden = o2ipa_num_from_map(summary_vals, "objective_burden"), objective_ploidy = o2ipa_num_from_map(summary_vals,
                "objective_ploidy"), burden_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw"),
            ploidy_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw"), runtime = o2ipa_num_from_map(summary_vals,
                "runtime"), parameter_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "best_params.tsv")))
                file.path(seed_dir, "best_params.tsv")
            else NA_character_, config_file = if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "fit_config.rds")))
                file.path(seed_dir, "fit_config.rds")
            else NA_character_, visualization_available = !is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "viz_status.log")), failure_reason = if (fit_success)
                NA_character_
            else paste(failure_parts, collapse = ";"), boundary_risk = boundary_risk, number_of_parameters_near_boundary = n_near,
            stringsAsFactors = FALSE)
    }
    manifest <- do.call(rbind, manifest_rows)
    params_long <- do.call(rbind, param_rows)
    list(manifest = manifest, params_long = params_long, boundary_long = boundary_long, summary_tab = summary_tab, matrix_tab = matrix_tab)
}

o2ipa_params_wide <- function(params_long, value_col = "value") {
    seeds <- unique(params_long$seed_id)
    params <- o2ipa_target_params()
    mat <- matrix(NA_real_, nrow = length(seeds), ncol = length(params), dimnames = list(seeds, params))
    for (i in seq_len(nrow(params_long))) {
        mat[params_long$seed_id[[i]], params_long$parameter[[i]]] <- params_long[[value_col]][[i]]
    }
    as.data.frame(mat, check.names = FALSE)
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

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
        y
    else x
}

o2pr_first_seed_cfg <- function(manifest) {
    for (p in manifest$config_file) {
        if (!is.na(p) && file.exists(p)) {
            cfg <- tryCatch(readRDS(p), error = function(e) NULL)
            if (!is.null(cfg))
                return(cfg)
        }
    }
    list(N_MIN = 22L, N_MAX = 154L, N_UNIT = 22L, DT = 0.05, start_with = "chr_number")
}

o2pr_cfg_metadata <- function(cfg) {
    data.frame(metric = c("Nmin", "Nmax", "N_UNIT", "DT", "start_with", "boundary_default", "o2_burden_feedback", "O2_growth",
        "ploidy_O2_death", "o2_S0_upper_bound", "o2_Nref", "trajectory_value_semantics"), value = c(cfg$N_MIN %||% NA, cfg$N_MAX %||%
        NA, cfg$N_UNIT %||% NA, cfg$DT %||% NA, cfg$start_with %||% NA, cfg$boundary %||% "drop", cfg$o2_burden_feedback %||%
        NA, cfg$O2_growth %||% NA, cfg$ploidy_O2_death %||% NA, cfg$o2_S0_upper_bound %||% NA, cfg$o2_Nref %||% NA, if (identical(as.character(cfg$start_with %||%
        "ploidy"), "chr_number")) {
        "weighted mean chromosome number N"
    } else {
        "weighted mean ploidy converted to N when needed"
    }), stringsAsFactors = FALSE)
}

o2pr_build_G <- function(model_env, cfg, run_params, O2) {
    fn <- get("cpp_o2simps_build_G_for_o2_triplet", envir = model_env, inherits = TRUE)
    tri <- fn(O2 = as.numeric(O2), O2_crit = as.numeric(run_params$O2_crit %||% cfg$o2_crit_init %||% 1), N0min = as.integer(cfg$N_MIN %||%
        22L), N0max = as.integer(cfg$N_MAX %||% 154L), N1min = as.integer(cfg$N_MIN %||% 22L), N1max = as.integer(cfg$N_MAX %||%
        154L), lam_max = as.numeric(run_params$lam_max), p_mis_base = as.numeric(run_params$p_mis_base %||% cfg$p_mis_base %||%
        1e-05), p_misseg = as.numeric(run_params$p_misseg %||% 0), k_o_mis = as.numeric(run_params$k_o_mis %||% 50), p_wgd = as.numeric(run_params$p_wgd %||%
        0), boundary = as.character(run_params$boundary %||% "drop"), eps_tail = 1e-08, buffer_smax = as.numeric(run_params$buffer_smax %||%
        1), buffer_beta = as.numeric(run_params$buffer_beta %||% 0), buffer_n_exp = as.numeric(run_params$buffer_n_exp %||%
        1), N_unit = as.integer(cfg$N_UNIT %||% 22L), beta_size = 0, O2_growth = isTRUE(run_params$O2_growth %||% cfg$O2_growth %||%
        TRUE), alpha_o2 = as.numeric(run_params$alpha_o2 %||% 0), gamma_growth = as.numeric(run_params$gamma_growth %||%
        1), mu_hp = as.numeric(run_params$mu_hp %||% 0), gamma_mu = as.numeric(run_params$gamma_mu %||% 1), n_O = as.numeric(run_params$n_O %||%
        1), ploidy_O2_death = as.character(run_params$ploidy_O2_death %||% cfg$ploidy_O2_death %||% "diploid_NULL"))
    G <- Matrix::sparseMatrix(i = as.integer(tri$i), j = as.integer(tri$j), x = as.numeric(tri$x), dims = c(as.integer(tri$nrow),
        as.integer(tri$ncol)), repr = "C")
    attr(G, "triplet") <- tri
    G
}

o2pr_run_params_from_vec <- function(vec, cfg) {
    rp <- as.list(vec)
    rp$o2_min <- if ("o2_min" %in% names(vec))
        vec[["o2_min"]]
    else (cfg$o2_min %||% 0)
    rp$o2_S0_upper_bound <- cfg$o2_S0_upper_bound %||% 5
    rp$o2_Nref <- cfg$o2_Nref %||% 1e+06
    rp$O2_growth <- cfg$O2_growth %||% TRUE
    rp$ploidy_O2_death <- cfg$ploidy_O2_death %||% "ploidy_related"
    rp$boundary <- cfg$boundary %||% "drop"
    rp
}
