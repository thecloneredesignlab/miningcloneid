#!/usr/bin/env Rscript

# Shared fixed-O2 fitting-parameter and model-loading helpers used by
# simulation/o2/fixed_o2/run_fixed_o2_simulation.R and analyses that source it.
suppressPackageStartupMessages({
  if (requireNamespace("Matrix", quietly = TRUE)) library(Matrix)
})

fixed_o2_canonical_utils_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_fixed_o2_utils.R"
  ]
  if (length(own)) {
    dirname(own[[length(own)]])
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})
source(
  file.path(
    fixed_o2_canonical_utils_dir,
    "o2_supply_demand_map_postfit_input_utils.R"
  ),
  local = TRUE
)
source(
  file.path(
    fixed_o2_canonical_utils_dir,
    "o2_supply_demand_map_postfit_probe_utils.R"
  ),
  local = TRUE
)
rm(fixed_o2_canonical_utils_dir)

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

o2ipa_path_is_absolute <- function(path) {
    grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

o2ipa_resolve_path <- function(path, base_dir = NULL, must_work = FALSE) {
    if (is.null(path) || !length(path))
        return(NA_character_)
    txt <- trimws(as.character(path[[1]]))
    if (!nzchar(txt) || is.na(txt) || txt %in% c("NA", "NaN", "NULL"))
        return(NA_character_)
    txt <- path.expand(txt)
    if (!o2ipa_path_is_absolute(txt) && !is.null(base_dir) && nzchar(base_dir)) {
        txt <- file.path(base_dir, txt)
    }
    normalizePath(txt, mustWork = must_work)
}

o2ipa_first_present_col <- function(tab, candidates) {
    hit <- candidates[candidates %in% names(tab)]
    if (length(hit))
        hit[[1]]
    else NA_character_
}

o2ipa_manifest_metric <- function(row, keys, default = NA_real_) {
    key <- o2ipa_first_present_col(row, keys)
    if (is.na(key))
        return(default)
    val <- suppressWarnings(as.numeric(row[[key]][[1]]))
    if (is.finite(val)) val else default
}

o2ipa_manifest_string <- function(row, keys, default = NA_character_) {
    key <- o2ipa_first_present_col(row, keys)
    if (is.na(key))
        return(default)
    val <- as.character(row[[key]][[1]])
    if (!is.na(val) && nzchar(val) && !val %in% c("NA", "NaN", "NULL")) val else default
}

o2ipa_read_seed_manifest <- function(seed_manifest) {
    manifest_path <- o2ipa_resolve_path(seed_manifest, must_work = TRUE)
    tab <- o2ipa_read_tsv(manifest_path)
    if (!"seed_id" %in% names(tab)) {
        if ("seed" %in% names(tab)) {
            tab$seed_id <- o2ipa_norm_seed(tab$seed)
        } else {
            stop("seed_manifest must contain seed_id or seed column: ", manifest_path)
        }
    }
    tab$seed_id <- o2ipa_norm_seed(tab$seed_id)
    if (anyDuplicated(tab$seed_id)) {
        stop("seed_manifest contains duplicated seed_id values: ", paste(unique(tab$seed_id[duplicated(tab$seed_id)]), collapse = ", "))
    }
    base_dir <- dirname(manifest_path)
    seed_dir_col <- o2ipa_first_present_col(tab, c("seed_dir", "source_seed_dir", "joint_seed_dir"))
    param_col <- o2ipa_first_present_col(tab, c("parameter_file", "best_params_file", "best_params_path"))
    summary_col <- o2ipa_first_present_col(tab, c("fit_summary_file", "summary_file", "fit_summary_path"))
    config_col <- o2ipa_first_present_col(tab, c("config_file", "fit_config_file", "fit_config_path"))
    if (!is.na(seed_dir_col)) {
        tab$seed_dir <- vapply(tab[[seed_dir_col]], o2ipa_resolve_path, character(1), base_dir = base_dir, must_work = FALSE)
    } else if (!"seed_dir" %in% names(tab)) {
        tab$seed_dir <- NA_character_
    }
    if (!is.na(param_col)) {
        tab$parameter_file <- vapply(tab[[param_col]], o2ipa_resolve_path, character(1), base_dir = base_dir, must_work = FALSE)
    } else if (!"parameter_file" %in% names(tab)) {
        tab$parameter_file <- ifelse(!is.na(tab$seed_dir) & nzchar(tab$seed_dir), file.path(tab$seed_dir, "best_params.tsv"), NA_character_)
    }
    if (!is.na(summary_col)) {
        tab$fit_summary_file <- vapply(tab[[summary_col]], o2ipa_resolve_path, character(1), base_dir = base_dir, must_work = FALSE)
    } else if (!"fit_summary_file" %in% names(tab)) {
        tab$fit_summary_file <- ifelse(!is.na(tab$seed_dir) & nzchar(tab$seed_dir), file.path(tab$seed_dir, "fit_summary.tsv"), NA_character_)
    }
    if (!is.na(config_col)) {
        tab$config_file <- vapply(tab[[config_col]], o2ipa_resolve_path, character(1), base_dir = base_dir, must_work = FALSE)
    } else if (!"config_file" %in% names(tab)) {
        tab$config_file <- ifelse(!is.na(tab$seed_dir) & nzchar(tab$seed_dir), file.path(tab$seed_dir, "fit_config.rds"), NA_character_)
    }
    attr(tab, "source_file") <- manifest_path
    tab
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










o2ipa_extract_all_params <- function(seed_id, seed_dir, summary_tab, matrix_tab, parameter_file = NA_character_) {
    best_path <- if (!is.na(parameter_file) && nzchar(parameter_file)) {
        parameter_file
    } else if (!is.na(seed_dir) && nzchar(seed_dir)) {
        file.path(seed_dir, "best_params.tsv")
    } else {
        NA_character_
    }
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


o2ipa_collect_seed_inputs <- function(run_dir, objective_source = "auto", seed_manifest = NULL) {
    manifest0 <- if (is.null(seed_manifest) || !length(seed_manifest) || is.na(seed_manifest[[1]]) || !nzchar(as.character(seed_manifest[[1]]))) {
        o2ipa_discover_seeds(run_dir)
    } else {
        o2ipa_read_seed_manifest(seed_manifest)
    }
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
        parameter_file <- if ("parameter_file" %in% names(manifest0))
            manifest0$parameter_file[[i]]
        else NA_character_
        if (is.na(parameter_file) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "best_params.tsv")))
            parameter_file <- file.path(seed_dir, "best_params.tsv")
        fit_summary_path <- if ("fit_summary_file" %in% names(manifest0) && !is.na(manifest0$fit_summary_file[[i]]) && nzchar(manifest0$fit_summary_file[[i]])) {
            manifest0$fit_summary_file[[i]]
        } else if (nzchar(seed_dir)) {
            file.path(seed_dir, "fit_summary.tsv")
        } else {
            NA_character_
        }
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
        manifest_row <- manifest0[i, , drop = FALSE]
        for (key in c("objective", "objective_total", "objective_data", "objective_burden", "objective_ploidy",
            "runtime", "objective_burden_neg2loglik_raw", "objective_ploidy_neg2loglik_raw")) {
            if (key %in% names(manifest_row) && !is.na(manifest_row[[key]][[1]]) && nzchar(as.character(manifest_row[[key]][[1]]))) {
                summary_vals[[key]] <- as.character(manifest_row[[key]][[1]])
            }
        }
        obj <- o2ipa_choose_objective(summary_vals, objective_source = objective_source)
        params_long <- o2ipa_extract_all_params(seed_id, seed_dir, summary_tab, matrix_tab, parameter_file = parameter_file)
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
        config_file <- if ("config_file" %in% names(manifest0) && !is.na(manifest0$config_file[[i]]) && nzchar(manifest0$config_file[[i]])) {
            manifest0$config_file[[i]]
        } else if (!is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir, "fit_config.rds"))) {
            file.path(seed_dir, "fit_config.rds")
        } else {
            NA_character_
        }
        manifest_rows[[i]] <- data.frame(seed_id = seed_id, seed_dir = if (nzchar(seed_dir))
            normalizePath(seed_dir, mustWork = FALSE)
        else NA_character_, fit_success = fit_success, convergence_status = convergence, objective = obj$value, objective_source = obj$source,
            objective_total = o2ipa_num_from_map(summary_vals, "objective_total"), objective_data = o2ipa_num_from_map(summary_vals,
                "objective_data"), objective_burden = o2ipa_num_from_map(summary_vals, "objective_burden"), objective_ploidy = o2ipa_num_from_map(summary_vals,
                "objective_ploidy"), burden_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_burden_neg2loglik_raw"),
            ploidy_neg2loglik_raw = o2ipa_num_from_map(summary_vals, "objective_ploidy_neg2loglik_raw"), runtime = o2ipa_num_from_map(summary_vals,
                "runtime"), parameter_file = if (!is.na(parameter_file) && nzchar(parameter_file) && file.exists(parameter_file))
                parameter_file
            else NA_character_, fit_summary_file = if (!is.na(fit_summary_path) && nzchar(fit_summary_path) && file.exists(fit_summary_path))
                fit_summary_path
            else NA_character_, config_file = if (!is.na(config_file) && nzchar(config_file) && file.exists(config_file))
                config_file
            else NA_character_, visualization_available = !is.na(seed_dir) && nzchar(seed_dir) && file.exists(file.path(seed_dir,
                "viz_status.log")), failure_reason = if (fit_success)
                NA_character_
            else paste(failure_parts, collapse = ";"), boundary_risk = boundary_risk, number_of_parameters_near_boundary = n_near,
            stringsAsFactors = FALSE)
        extra_cols <- setdiff(names(manifest_row), names(manifest_rows[[i]]))
        for (col in extra_cols) manifest_rows[[i]][[col]] <- manifest_row[[col]][[1]]
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



`%||%` <- function(x, y) {
    if (is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
        y
    else x
}





bpf_default_fixed_o2_grid <-
function() {
    c(0, 0.1, 0.5, 1, 2, 5)
}

bpf_default_dense_attractor_o2_grid <-
function() {
    seq(0, 5, by = 0.05)
}

bpf_o2_key <-
function(x) {
    vapply(x, function(xx) format(signif(as.numeric(xx), 12), scientific = FALSE, trim = TRUE), character(1))
}

bpf_o2_slug <-
function(x) {
    key <- bpf_o2_key(x)
    key <- gsub("-", "minus", key, fixed = TRUE)
    key <- gsub("[^0-9A-Za-z]+", "p", key)
    key <- gsub("^p+|p+$", "", key)
    ifelse(nzchar(key), key, "NA")
}

fixed_o2_shared_utils_dir <-
function() {
    frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
        ofile <- env$ofile
        if (is.null(ofile))
            ""
        else normalizePath(ofile, mustWork = FALSE)
    }, character(1)))
    if (length(frame_files)) {
        return(dirname(frame_files[[length(frame_files)]]))
    }
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)))
    }
    normalizePath(getwd(), mustWork = FALSE)
}
