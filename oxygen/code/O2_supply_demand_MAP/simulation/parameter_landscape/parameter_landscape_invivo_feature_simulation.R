#!/usr/bin/env Rscript

# Model-backed in-vivo feature reconstruction for the parameter-landscape simulation layer.
.o2pl_feature_dir <- local({
  frames <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile; if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frames[basename(frames) == "parameter_landscape_invivo_feature_simulation.R"]
  if (length(own)) dirname(own[[length(own)]]) else normalizePath(getwd(), mustWork = FALSE)
})
source(file.path(.o2pl_feature_dir, "parameter_landscape_simulation_utils.R"), local = environment(), chdir = TRUE)

default_invivo_data_dir <- function() {
    oxygen_dir <- default_oxygen_dir()
    candidates <- c(file.path(oxygen_dir, "..", "data", "InVivoData_Gemcitabine"), file.path(oxygen_dir, "data", "InVivoData_Gemcitabine"),
        file.path(dirname(oxygen_dir), "data", "InVivoData_Gemcitabine"))
    for (candidate in candidates) {
        if (file.exists(file.path(candidate, "dt_Gem_VT_20260209_v5.xlsx")) && file.exists(file.path(candidate, "all_ploidy.csv"))) {
            return(normalizePath(candidate, mustWork = FALSE))
        }
    }
    normalizePath(candidates[[1L]], mustWork = FALSE)
}

invivo_data_paths <- function(data_dir = default_invivo_data_dir()) {
    data_dir <- normalizePath(path.expand(data_dir), mustWork = FALSE)
    dt_path <- file.path(data_dir, "dt_Gem_VT_20260209_v5.xlsx")
    ploidy_path <- file.path(data_dir, "all_ploidy.csv")
    if (!file.exists(dt_path))
        stop("Missing in vivo tumor-burden file: ", dt_path)
    if (!file.exists(ploidy_path))
        stop("Missing in vivo ploidy file: ", ploidy_path)
    list(dt_path = dt_path, ploidy_path = ploidy_path)
}

invivo_simulation_env <- local({
    cache <- NULL
    function() {
        if (!is.null(cache))
            return(cache)
        workflow_root <- normalizePath(file.path(default_oxygen_dir(), "code", "O2_supply_demand_MAP"), mustWork = FALSE)
        shared_path <- file.path(workflow_root, "util", "o2_supply_demand_map_shared.R")
        common_path <- file.path(workflow_root, "util", "o2_supply_demand_map_common_semantics.R")
        model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")
        simulation_utils_path <- file.path(workflow_root, "simulation", "invivo", "o2_supply_demand_map_invivo_simulation_utils.R")
        if (!file.exists(shared_path))
            stop("Shared utility script not found: ", shared_path)
        if (!file.exists(common_path))
            stop("Common semantics script not found: ", common_path)
        if (!file.exists(model_path))
            stop("Model script not found: ", model_path)
        if (!file.exists(simulation_utils_path)) {
            stop("In vivo simulation utilities not found: ", simulation_utils_path)
        }
        env <- new.env(parent = globalenv())
        env$commandArgs <- function(trailingOnly = FALSE) character(0)
        sys.source(shared_path, envir = globalenv(), chdir = TRUE)
        sys.source(common_path, envir = globalenv(), chdir = TRUE)
        sys.source(shared_path, envir = env, chdir = TRUE)
        sys.source(common_path, envir = env, chdir = TRUE)
        Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = dirname(model_path))
        sys.source(model_path, envir = env, chdir = TRUE)
        sys.source(simulation_utils_path, envir = env, chdir = TRUE)
        required <- c("normalize_cfg_for_simulation", "read_run_params", "prepare_data", "simulate_one_full_horizon", ".lambda_eff_of_O2",
            ".mu_eff_of_O2", ".pmisseg_of_O2", ".p_wgd_of_O2", ".pr_delta_vec")
        missing <- required[!vapply(required, exists, logical(1), envir = env, inherits = TRUE)]
        if (length(missing))
            stop("In vivo simulation environment is missing: ", paste(missing, collapse = ", "))
        cache <<- env
        cache
    }
})

compute_misseg_death_rate <- function(lambda, O2, N, run_params, cfg = NULL, model_env = invivo_simulation_env()) {
    lambda <- suppressWarnings(as.numeric(lambda))
    O2 <- suppressWarnings(as.numeric(O2))
    N <- suppressWarnings(as.numeric(N))
    if (!(length(lambda) == length(O2) && length(O2) == length(N))) {
        stop("lambda, O2, and N must have equal lengths.")
    }
    p_mis <- get(".pmisseg_of_O2", envir = model_env, inherits = TRUE)(O2 = O2, run_params = run_params, N = N)
    p_wgd <- get(".p_wgd_of_O2", envir = model_env, inherits = TRUE)(O2 = O2, run_params = run_params)
    pr_delta <- get(".pr_delta_vec", envir = model_env, inherits = TRUE)
    buffer_smax <- as.numeric(run_params$buffer_smax %||% 1)
    if (!is.finite(buffer_smax))
        buffer_smax <- 1
    buffer_beta <- as.numeric(run_params$buffer_beta %||% 0)
    if (!is.finite(buffer_beta))
        buffer_beta <- 0
    buffer_n_exp <- as.numeric(run_params$buffer_n_exp %||% 1)
    if (!is.finite(buffer_n_exp))
        buffer_n_exp <- 1
    n_unit <- as.integer((cfg$N_UNIT %||% cfg$N_unit %||% run_params$N_UNIT %||% 22L)[[1L]])
    if (!is.finite(n_unit) || n_unit <= 0L)
        n_unit <- 22L
    cache <- new.env(parent = emptyenv())
    mass_dropped <- vapply(seq_along(N), function(i) {
        if (!is.finite(N[[i]]) || !is.finite(p_mis[[i]]))
            return(NA_real_)
        n_i <- as.integer(round(N[[i]]))
        p_i <- pmin(pmax(as.numeric(p_mis[[i]]), 0), 1)
        key <- paste(n_i, signif(p_i, 12), sep = "|")
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }
        delta <- pr_delta(N = n_i, p = p_i, eps_tail = 1e-08, N_unit = n_unit, buffer_smax = buffer_smax, buffer_beta = buffer_beta,
            buffer_n_exp = buffer_n_exp)
        value <- as.numeric(attr(delta, "mass_dropped"))
        if (!is.finite(value))
            value <- NA_real_
        assign(key, value, envir = cache)
        value
    }, numeric(1))
    dead_daughters_per_division <- pmin(pmax(2 * mass_dropped, 0), 2)
    p_wgd <- pmin(pmax(as.numeric(p_wgd), 0), 1)
    lambda * (1 - p_wgd) * dead_daughters_per_division
}

time_weighted_mean <- function(day, value) {
    day <- suppressWarnings(as.numeric(day))
    value <- suppressWarnings(as.numeric(value))
    keep <- is.finite(day) & is.finite(value)
    day <- day[keep]
    value <- value[keep]
    if (!length(day))
        return(NA_real_)
    ord <- order(day)
    day <- day[ord]
    value <- value[ord]
    day_unique <- sort(unique(day))
    if (length(day_unique) < length(day)) {
        value <- vapply(day_unique, function(d) mean(value[day == d], na.rm = TRUE), numeric(1))
        day <- day_unique
    }
    if (length(day) == 1L || max(day) <= min(day))
        return(mean(value, na.rm = TRUE))
    widths <- diff(day)
    integral <- sum(widths * (head(value, -1L) + tail(value, -1L))/2, na.rm = TRUE)
    integral/(max(day) - min(day))
}

compute_invivo_rate_summary_for_sim <- function(sim, run_params, cfg, model_env, horizon_day = 100) {
    burden <- sim$burden
    ploidy <- sim$ploidy
    if (is.null(burden) || is.null(ploidy) || !nrow(burden) || !nrow(ploidy)) {
        return(NULL)
    }
    burden <- burden[is.finite(burden$day) & burden$day <= horizon_day + 1e-09, , drop = FALSE]
    ploidy <- ploidy[is.finite(ploidy$day) & ploidy$day <= horizon_day + 1e-09, , drop = FALSE]
    if (!nrow(burden) || !nrow(ploidy))
        return(NULL)
    mean_o2 <- time_weighted_mean(burden$day, burden$pred_o2_pct)
    o2_by_day <- stats::setNames(as.numeric(burden$pred_o2_pct), as.character(as.numeric(burden$day)))
    ploidy$o2_pct <- as.numeric(o2_by_day[as.character(as.numeric(ploidy$day))])
    keep <- is.finite(ploidy$o2_pct) & is.finite(ploidy$N) & is.finite(ploidy$fraction)
    ploidy <- ploidy[keep, , drop = FALSE]
    if (!nrow(ploidy))
        return(NULL)
    lambda <- get(".lambda_eff_of_O2", envir = model_env, inherits = TRUE)(O2 = ploidy$o2_pct, run_params = run_params, N = ploidy$N,
        O2_growth = isTRUE(cfg$O2_growth %||% TRUE))
    mu <- get(".mu_eff_of_O2", envir = model_env, inherits = TRUE)(O2 = ploidy$o2_pct, run_params = run_params, N = ploidy$N)
    misseg_death <- compute_misseg_death_rate(lambda = lambda, O2 = ploidy$o2_pct, N = ploidy$N, run_params = run_params,
        cfg = cfg, model_env = model_env)
    rate_df <- data.frame(day = as.numeric(ploidy$day), fraction = pmax(as.numeric(ploidy$fraction), 0), lambda = as.numeric(lambda),
        mu = as.numeric(mu), misseg_death = as.numeric(misseg_death), stringsAsFactors = FALSE)
    rate_df <- rate_df[is.finite(rate_df$day) & is.finite(rate_df$fraction) & is.finite(rate_df$lambda) & is.finite(rate_df$mu) &
        is.finite(rate_df$misseg_death), , drop = FALSE]
    if (!nrow(rate_df))
        return(NULL)
    by_day <- split(rate_df, rate_df$day)
    day_rates <- do.call(rbind, lapply(by_day, function(df) {
        w <- pmax(df$fraction, 0)
        w_sum <- sum(w, na.rm = TRUE)
        if (!is.finite(w_sum) || w_sum <= 0)
            return(NULL)
        lambda_mean <- sum(w * df$lambda, na.rm = TRUE)/w_sum
        mu_mean <- sum(w * df$mu, na.rm = TRUE)/w_sum
        misseg_death_mean <- sum(w * df$misseg_death, na.rm = TRUE)/w_sum
        data.frame(day = as.numeric(df$day[[1L]]), net_growth = lambda_mean - mu_mean, turnover_rate = lambda_mean + mu_mean,
            misseg_death_rate = misseg_death_mean, net_growth_with_misseg_death = lambda_mean - mu_mean - misseg_death_mean,
            turnover_rate_with_misseg_death = lambda_mean + mu_mean + misseg_death_mean, stringsAsFactors = FALSE)
    }))
    if (is.null(day_rates) || !nrow(day_rates))
        return(NULL)
    data.frame(mean_net_growth_0_100d = time_weighted_mean(day_rates$day, day_rates$net_growth), mean_turnover_rate_0_100d = time_weighted_mean(day_rates$day,
        day_rates$turnover_rate), mean_O2_0_100d = mean_o2, mean_misseg_death_rate_0_100d = time_weighted_mean(day_rates$day,
        day_rates$misseg_death_rate), mean_net_growth_with_misseg_death_0_100d = time_weighted_mean(day_rates$day, day_rates$net_growth_with_misseg_death),
        mean_turnover_rate_with_misseg_death_0_100d = time_weighted_mean(day_rates$day, day_rates$turnover_rate_with_misseg_death),
        stringsAsFactors = FALSE)
}

select_invivo_representative_scenarios <- function(scenarios) {
    cohorts <- vapply(scenarios, function(sc) as.character(sc$cohort %||% ""), character(1))
    reps <- lapply(c("2N", "4N"), function(cohort_id) {
        idx <- which(cohorts == cohort_id)
        if (!length(idx))
            stop("No ", cohort_id, " scenario found for growth/turnover calculation.")
        scenarios[[idx[[1L]]]]
    })
    names(reps) <- c("2N", "4N")
    reps
}

compute_invivo_growth_turnover_one_seed <- function(seed_dir, data_paths, model_env = invivo_simulation_env(), horizon_day = 100,
    report_dt = 1) {
    seed <- seed_from_dir(seed_dir)
    cfg_path <- file.path(seed_dir, "fit_config.rds")
    if (!file.exists(cfg_path))
        stop("Missing fit_config.rds: ", cfg_path)
    cfg <- readRDS(cfg_path)
    cfg <- get("normalize_cfg_for_simulation", envir = model_env, inherits = TRUE)(cfg)
    run_params <- get("read_run_params", envir = model_env, inherits = TRUE)(seed_dir, cfg = cfg)
    scenarios <- get("prepare_data", envir = model_env, inherits = TRUE)(data_paths$dt_path, data_paths$ploidy_path, cfg)
    scenarios <- select_invivo_representative_scenarios(scenarios)
    sim_fun <- get("simulate_one_full_horizon", envir = model_env, inherits = TRUE)
    scenario_rows <- vector("list", length(scenarios))
    for (i in seq_along(scenarios)) {
        sc <- scenarios[[i]]
        sim <- sim_fun(run_params = run_params, scenario = sc, cfg = cfg, horizon_day = horizon_day, report_dt = report_dt)
        rates <- compute_invivo_rate_summary_for_sim(sim = sim, run_params = run_params, cfg = cfg, model_env = model_env,
            horizon_day = horizon_day)
        if (is.null(rates))
            next
        scenario_rows[[i]] <- data.frame(seed = seed, cohort = as.character(sc$cohort), harvest = as.character(sc$harvest),
            dose = as.numeric(sc$dose), rates, stringsAsFactors = FALSE)
    }
    scenario_df <- do.call(rbind, scenario_rows)
    if (is.null(scenario_df) || !nrow(scenario_df)) {
        stop("No rate summaries were generated for seed", seed)
    }
    cohort_levels <- c("2N", "4N")
    cohort_summary <- do.call(rbind, lapply(cohort_levels, function(cohort_id) {
        df <- scenario_df[scenario_df$cohort == cohort_id, , drop = FALSE]
        if (!nrow(df))
            return(NULL)
        data.frame(cohort = cohort_id, mean_net_growth_0_100d = mean(df$mean_net_growth_0_100d, na.rm = TRUE), mean_turnover_rate_0_100d = mean(df$mean_turnover_rate_0_100d,
            na.rm = TRUE), mean_O2_0_100d = mean(df$mean_O2_0_100d, na.rm = TRUE), mean_misseg_death_rate_0_100d = mean(df$mean_misseg_death_rate_0_100d,
            na.rm = TRUE), mean_net_growth_with_misseg_death_0_100d = mean(df$mean_net_growth_with_misseg_death_0_100d, na.rm = TRUE),
            mean_turnover_rate_with_misseg_death_0_100d = mean(df$mean_turnover_rate_with_misseg_death_0_100d, na.rm = TRUE),
            n_scenarios = nrow(df), stringsAsFactors = FALSE)
    }))
    if (is.null(cohort_summary) || !nrow(cohort_summary)) {
        stop("No 2N/4N cohort summaries were generated for seed", seed)
    }
    value_for <- function(cohort_id, col) {
        hit <- cohort_summary[cohort_summary$cohort == cohort_id, col, drop = TRUE]
        if (length(hit))
            as.numeric(hit[[1L]])
        else NA_real_
    }
    n_for <- function(cohort_id) {
        hit <- cohort_summary[cohort_summary$cohort == cohort_id, "n_scenarios", drop = TRUE]
        if (length(hit))
            as.integer(hit[[1L]])
        else 0L
    }
    net_2n <- value_for("2N", "mean_net_growth_0_100d")
    net_4n <- value_for("4N", "mean_net_growth_0_100d")
    turnover_2n <- value_for("2N", "mean_turnover_rate_0_100d")
    turnover_4n <- value_for("4N", "mean_turnover_rate_0_100d")
    o2_2n <- value_for("2N", "mean_O2_0_100d")
    o2_4n <- value_for("4N", "mean_O2_0_100d")
    misseg_death_2n <- value_for("2N", "mean_misseg_death_rate_0_100d")
    misseg_death_4n <- value_for("4N", "mean_misseg_death_rate_0_100d")
    net_with_misseg_death_2n <- value_for("2N", "mean_net_growth_with_misseg_death_0_100d")
    net_with_misseg_death_4n <- value_for("4N", "mean_net_growth_with_misseg_death_0_100d")
    turnover_with_misseg_death_2n <- value_for("2N", "mean_turnover_rate_with_misseg_death_0_100d")
    turnover_with_misseg_death_4n <- value_for("4N", "mean_turnover_rate_with_misseg_death_0_100d")
    data.frame(seed = seed, mean_net_growth_0_100d_2N = net_2n, mean_net_growth_0_100d_4N = net_4n, mean_net_growth_0_100d = mean(c(net_2n,
        net_4n), na.rm = TRUE), mean_turnover_rate_0_100d_2N = turnover_2n, mean_turnover_rate_0_100d_4N = turnover_4n, mean_turnover_rate_0_100d = mean(c(turnover_2n,
        turnover_4n), na.rm = TRUE), mean_O2_0_100d_2N = o2_2n, mean_O2_0_100d_4N = o2_4n, mean_O2_0_100d = mean(c(o2_2n,
        o2_4n), na.rm = TRUE), mean_misseg_death_rate_0_100d_2N = misseg_death_2n, mean_misseg_death_rate_0_100d_4N = misseg_death_4n,
        mean_misseg_death_rate_0_100d = mean(c(misseg_death_2n, misseg_death_4n), na.rm = TRUE), mean_net_growth_with_misseg_death_0_100d_2N = net_with_misseg_death_2n,
        mean_net_growth_with_misseg_death_0_100d_4N = net_with_misseg_death_4n, mean_net_growth_with_misseg_death_0_100d = mean(c(net_with_misseg_death_2n,
            net_with_misseg_death_4n), na.rm = TRUE), mean_turnover_rate_with_misseg_death_0_100d_2N = turnover_with_misseg_death_2n,
        mean_turnover_rate_with_misseg_death_0_100d_4N = turnover_with_misseg_death_4n, mean_turnover_rate_with_misseg_death_0_100d = mean(c(turnover_with_misseg_death_2n,
            turnover_with_misseg_death_4n), na.rm = TRUE), n_rate_scenarios_2N = n_for("2N"), n_rate_scenarios_4N = n_for("4N"),
        n_rate_scenarios = nrow(scenario_df), rate_horizon_day = as.numeric(horizon_day), rate_report_dt = as.numeric(report_dt),
        stringsAsFactors = FALSE)
}

paper_generate_invivo_growth_turnover_table <- function(input_dir = default_dataset_input_dir("invivo"), best_csv = paper_best_params_csv("invivo"),
    output_csv = file.path(paper_tables_dir("invivo"), "invivo_best_params_growth_turnover_100d.csv"), data_dir = default_invivo_data_dir(),
    max_seeds = NA_integer_, horizon_day = 100, report_dt = 1) {
    input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
    best_csv <- normalizePath(path.expand(best_csv), mustWork = FALSE)
    output_csv <- normalizePath(path.expand(output_csv), mustWork = FALSE)
    max_seeds <- as_int(max_seeds, NA_integer_)
    horizon_day <- as_num(horizon_day, 100)
    report_dt <- as_num(report_dt, 1)
    if (!dir.exists(input_dir))
        stop("Input directory does not exist: ", input_dir)
    if (!file.exists(best_csv))
        stop("Missing best-parameter CSV: ", best_csv)
    if (!is.finite(horizon_day) || horizon_day <= 0)
        stop("horizon_day must be > 0.")
    if (!is.finite(report_dt) || report_dt <= 0)
        stop("report_dt must be > 0.")
    dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
    best_df <- read_csv_plain(best_csv)
    if (!"seed" %in% names(best_df))
        stop("Best-parameter CSV must contain a seed column.")
    best_df$seed <- as.integer(best_df$seed)
    seed_dirs <- list_seed_dirs(input_dir)
    if (!length(seed_dirs))
        stop("No seed directories found under: ", input_dir)
    if (is.finite(max_seeds) && !is.na(max_seeds) && max_seeds > 0L) {
        seed_dirs <- seed_dirs[seq_len(min(length(seed_dirs), max_seeds))]
    }
    seed_ids <- vapply(seed_dirs, seed_from_dir, integer(1))
    missing_best <- setdiff(seed_ids, best_df$seed)
    if (length(missing_best)) {
        stop("Best-parameter CSV is missing seeds: ", paste(head(missing_best, 20), collapse = ", "))
    }
    data_paths <- invivo_data_paths(data_dir)
    model_env <- invivo_simulation_env()
    message("Computing in vivo 0-", horizon_day, " day growth/turnover metrics for ", length(seed_dirs), " seeds.")
    rows <- vector("list", length(seed_dirs))
    for (i in seq_along(seed_dirs)) {
        rows[[i]] <- compute_invivo_growth_turnover_one_seed(seed_dir = seed_dirs[[i]], data_paths = data_paths, model_env = model_env,
            horizon_day = horizon_day, report_dt = report_dt)
        if (i%%10L == 0L || i == length(seed_dirs)) {
            message("Computed growth/turnover metrics for ", i, "/", length(seed_dirs), " seeds.")
        }
    }
    metrics_df <- do.call(rbind, rows)
    best_subset <- best_df[match(metrics_df$seed, best_df$seed), , drop = FALSE]
    metric_cols <- setdiff(names(metrics_df), "seed")
    out <- cbind(best_subset, metrics_df[, metric_cols, drop = FALSE])
    out <- append_invivo_pred1000_ploidy_ratios(out, input_dir = input_dir)
    required_rate_cols <- c("mean_net_growth_0_100d", "mean_turnover_rate_0_100d", "mean_O2_0_100d", "mean_misseg_death_rate_0_100d",
        "mean_net_growth_with_misseg_death_0_100d", "mean_turnover_rate_with_misseg_death_0_100d")
    if (any(vapply(required_rate_cols, function(col) any(!is.finite(out[[col]])), logical(1)))) {
        stop("Computed growth/turnover table contains non-finite aggregate rate metrics.")
    }
    utils::write.csv(out, output_csv, quote = FALSE, row.names = FALSE)
    message("Wrote in vivo growth/turnover table: ", output_csv)
    invisible(output_csv)
}

append_invivo_pred1000_ploidy_ratios <- function(df, input_dir = default_dataset_input_dir("invivo"), target_day = 1000) {
    if (!"seed" %in% names(df))
        stop("Input table must contain a seed column.")
    input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
    seed_ids <- as.integer(df$seed)
    rows <- lapply(seed_ids, function(seed) {
        seed_dir <- file.path(input_dir, paste0("seed", seed))
        read_pred1000_ploidy_ratios(seed_dir, seed, target_day = target_day)
    })
    ratio_df <- do.call(rbind, rows)
    for (col in names(ratio_df)) {
        df[[col]] <- ratio_df[[col]]
    }
    ratio_cols <- names(ratio_df)
    if (any(vapply(ratio_cols, function(col) any(!is.finite(df[[col]])), logical(1)))) {
        stop("Computed pred1000 ploidy ratio columns contain non-finite values.")
    }
    df
}

paper_update_invivo_growth_turnover_ploidy_ratios <- function(input_dir = default_dataset_input_dir("invivo"), growth_turnover_csv = file.path(paper_tables_dir("invivo"),
    "invivo_best_params_growth_turnover_100d.csv"), target_day = 1000) {
    input_dir <- normalizePath(path.expand(input_dir), mustWork = FALSE)
    growth_turnover_csv <- normalizePath(path.expand(growth_turnover_csv), mustWork = FALSE)
    if (!dir.exists(input_dir))
        stop("Input directory does not exist: ", input_dir)
    if (!file.exists(growth_turnover_csv))
        stop("Missing growth/turnover CSV: ", growth_turnover_csv)
    df <- read_csv_plain(growth_turnover_csv)
    out <- append_invivo_pred1000_ploidy_ratios(df, input_dir = input_dir, target_day = target_day)
    utils::write.csv(out, growth_turnover_csv, quote = FALSE, row.names = FALSE)
    message("Updated in vivo growth/turnover table with pred1000 ploidy ratios: ", growth_turnover_csv)
    invisible(growth_turnover_csv)
}

rm(.o2pl_feature_dir)
