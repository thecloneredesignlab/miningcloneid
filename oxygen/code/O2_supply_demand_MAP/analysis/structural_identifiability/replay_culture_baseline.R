#!/usr/bin/env Rscript
# WP0 only: one forward evaluation from archived, full-precision parameters.
# Deliberately does not source a fit backend or call an optimizer.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: replay_culture_baseline.R <isolated-study-root>")
root <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(root, "replay")
if (file.exists(file.path(out, "replay_status.json"))) stop("Replay already completed; use a new study root.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
settings <- jsonlite::read_json(file.path(root, "manifests", "replay_settings.json"), simplifyVector = TRUE)
atol <- settings$absolute_tolerance
rtol <- settings$relative_tolerance
stopifnot(identical(settings$optimization, FALSE), identical(settings$passage_implementation, "v2"))
snapshot <- file.path(root, "source", "code", "oxygen", "code", "O2_supply_demand_MAP")
inputs <- file.path(root, "source", "inputs", "oxygen")
fit <- readRDS(file.path(root, "source", "seed144", "fit_result.rds"))
reference <- fit$best_components
cfg <- fit$cfg
original_parameter_path <- cfg$parameter_table
cfg$parameter_table <- file.path(root, "source", "seed144", "parameter_table_input.csv")
stopifnot(identical(reference$death_weight, 0), reference$n_death_observations == 0L)
Sys.setenv(MININGCLONEID_OXYGEN_CODE_DIR = file.path(snapshot, "model"),
           MININGCLONEID_RCPP_REBUILD = "TRUE")
started <- Sys.time()
source(file.path(snapshot, "model", "model_O2_supply_demand_MAP.R"), local = .GlobalEnv)
source(file.path(snapshot, "util", "o2_supply_demand_map_invitro_utils.R"), local = .GlobalEnv)
stopifnot(identical(INVITRO_PASSAGE_IMPLEMENTATION, "v2"))
objects <- ivt_load_fit_objects(
  repo_root = inputs,
  fit_objects_dir = file.path(inputs, "ploidyOxygen", "data", "fit_objects"),
  flow_csv_path = file.path(inputs, "data", "g0g1_ploidy_density_grid.csv"))
stopifnot(identical(objects$death_enabled, FALSE))
forward_started <- Sys.time()
result <- ivt_objective_components(
  run_params = fit$best_params, fit_objects = objects, cfg = cfg,
  fallback_max_passage_days = 14, growth_weight = 1, ploidy_weight = 1,
  flow_weight = 1, ploidy_prob_floor = 1e-12, flow_density_floor = 1e-12)
forward_seconds <- as.numeric(difftime(Sys.time(), forward_started, units = "secs"))
saveRDS(result, file.path(out, "replay_components.rds"))

checks <- list()
record <- function(scope, item, expected, actual, exact = FALSE) {
  shape_ok <- identical(dim(expected), dim(actual)) && length(expected) == length(actual)
  numeric <- is.numeric(expected) && is.numeric(actual)
  max_abs <- max_scaled <- NA_real_
  ok <- shape_ok
  if (ok && numeric) {
    ok <- identical(is.na(expected), is.na(actual)) && identical(is.finite(expected), is.finite(actual))
    special <- !is.na(expected) & !is.finite(expected)
    if (ok && any(special)) ok <- identical(expected[special], actual[special])
    keep <- is.finite(expected) & is.finite(actual)
    if (any(keep)) {
      diff <- abs(expected[keep] - actual[keep])
      max_abs <- max(diff)
      scale <- if (exact) rep(0, length(diff)) else atol + rtol * pmax(abs(expected[keep]), abs(actual[keep]))
      max_scaled <- if (exact) if (all(diff == 0)) 0 else Inf else max(diff / scale)
      ok <- ok && all(if (exact) diff == 0 else diff <= scale)
    }
  } else if (ok) ok <- identical(as.character(expected), as.character(actual))
  checks[[length(checks) + 1L]] <<- data.frame(
    scope = scope, item = item, n_expected = length(expected), n_actual = length(actual),
    max_absolute_difference = max_abs, max_scaled_difference = max_scaled,
    exact = exact, pass = ok, stringsAsFactors = FALSE)
}
compare_frame <- function(scope, a, b) {
  record(scope, "columns", names(a), names(b), TRUE)
  record(scope, "rows", as.integer(nrow(a)), as.integer(nrow(b)), TRUE)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) return(invisible(NULL))
  for (n in names(a)) record(scope, n, a[[n]], b[[n]], exact = n %in% c("selected_day", "passage_index", "grid_index"))
}
metrics <- c("objective", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik",
             "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum",
             "sigma_growth", "sigma_kary", "sigma_flow_ploidy", "n_growth", "n_growth_observed",
             "n_growth_missing_pred", "n_growth_negative_pred", "n_ploidy_passages", "n_kary_cells",
             "n_flow_passages", "n_flow_samples", "n_death_observations", "death_weight")
baseline <- do.call(rbind, lapply(metrics, function(n) {
  record("objective", n, reference[[n]], result[[n]], grepl("^n_|weight$", n))
  data.frame(metric = n, archived = reference[[n]], replay = result[[n]], difference = result[[n]] - reference[[n]])
}))
write_tsv(baseline, file.path(root, "audit", "baseline_replay.tsv"))
for (n in c("summary", "growth_df", "ploidy_df", "flow_df", "flow_overlay_df")) {
  compare_frame(n, reference[[n]], result[[n]])
  write_tsv(result[[n]], file.path(out, paste0(n, ".tsv")))
}

# A separate implementation of scoring arithmetic, using observations and
# selected predicted distributions, not returned component scores.
independent <- list(growth = numeric(), karyotype = numeric(), flow = numeric())
observations <- boundaries <- selections <- distributions <- per_sample <- list()
for (cohort in c("2N", "4N")) {
  run <- result[[paste0("run_", cohort)]]
  old <- reference[[paste0("run_", cohort)]]
  record(cohort, "segment_names", names(old$segment_results), names(run$segment_results), TRUE)
  record(cohort, "chromosome_grid", old$grid_pre, run$grid_pre, TRUE)
  for (i in seq_along(run$segment_results)) {
    s <- run$segment_results[[i]]
    a <- old$segment_results[[i]]
    seg <- s$segment
    scope <- paste(cohort, seg$segment_id, sep = "::")
    record(scope, "data_ids", a$segment$data_ids, seg$data_ids, TRUE)
    record(scope, "parent_segment_id", a$segment$parent_segment_id, seg$parent_segment_id, TRUE)
    record(scope, "observation_times", a$segment$obs_days_local, seg$obs_days_local, TRUE)
    for (n in names(s$selection)) record(scope, paste0("selection/", n), a$selection[[n]], s$selection[[n]], n %in% c("selected_index", "selected_day"))
    record(scope, "parent_boundary_present", !is.null(a$parent_boundary_selection), !is.null(s$parent_boundary_selection), TRUE)
    if (!is.null(s$parent_boundary_selection)) for (n in names(s$parent_boundary_selection)) {
      record(scope, paste0("propagation/", n), a$parent_boundary_selection[[n]], s$parent_boundary_selection[[n]], n %in% c("selected_index", "selected_day"))
    }
    for (n in setdiff(names(s$sim), c("cache_g_build", "cache_g_hit", "cache_g_hysteresis", "cache_bin_pct", "cache_hysteresis_pct"))) {
      record(scope, paste0("trajectory/", n), a$sim[[n]], s$sim[[n]])
    }
    candidate <- which(is.finite(seg$obs_days_local) & seg$obs_days_local > 0)
    if (!length(candidate)) candidate <- seq_along(seg$obs_days_local)
    target <- s$selection$target_live_cells
    delta <- if (is.finite(target) && target > 0) abs(s$sim$Ntot_live_obs[candidate] - target) else NA_real_
    selections[[length(selections) + 1L]] <- data.frame(
      cohort = cohort, segment_id = seg$segment_id, parent_segment_id = seg$parent_segment_id,
      oxygen_pct = seg$oxygen_pct, nominal_duration_days = seg$duration_days,
      selected_index = s$selection$selected_index, selected_day = s$selection$selected_day,
      target_live_cells = target, selected_live_cells = s$selection$selected_live_cells,
      candidate_days = paste(seg$obs_days_local[candidate], collapse = ","),
      nearest_count_ties = if (all(is.na(delta))) NA_integer_ else sum(delta == min(delta)),
      selection_rule = if (all(is.na(delta))) "last_candidate_no_positive_target" else "nearest_target_count_then_first_index")
    pb <- s$parent_boundary_selection
    if (!is.null(pb)) {
      ps <- run$segment_results[[seg$parent_segment_id]]
      boundaries[[length(boundaries) + 1L]] <- data.frame(
        cohort = cohort, child_segment_id = seg$segment_id, parent_segment_id = seg$parent_segment_id,
        parent_observation_day = ps$selection$selected_day,
        propagation_day = pb$selected_day, required_child_inoculum = pb$required_cells,
        available_parent_cells = pb$available_cells, supply_ratio = pb$supply_ratio,
        boundary_scale = pb$boundary_scale, reseed_mode = pb$reseed_mode,
        passage_executed = pb$passage_executed,
        differs_from_observation_day = pb$selected_day != ps$selection$selected_day)
    }
    p <- pmax(as.numeric(s$selection$selected_frac), 0)
    stopifnot(all(is.finite(p)), sum(p) > 0)
    p <- p / sum(p)
    distributions[[length(distributions) + 1L]] <- data.frame(
      cohort = cohort, segment_id = seg$segment_id, chromosome_number = run$grid_pre,
      selected_day = s$selection$selected_day, live_fraction = p)
    kernel <- outer(run$grid_pre, run$grid_pre,
                    function(observed, latent) stats::dnorm(observed, latent, fit$best_params$sigma_kary))
    cs <- colSums(kernel)
    cs[!is.finite(cs) | cs <= 0] <- 1
    smooth <- as.numeric(sweep(kernel, 2L, cs, "/") %*% p)
    names(smooth) <- as.character(run$grid_pre)
    for (pid in seg$data_ids) {
      obs <- ivt_observed_passage_summary(objects$fit_data[[pid]])
      fg <- obs$observed_flow
      present_flow <- !is.null(fg) && nrow(fg) > 0L
      observations[[length(observations) + 1L]] <- data.frame(
        observation_id = pid, cohort = cohort, segment_id = seg$segment_id,
        parent_segment_id = seg$parent_segment_id, passage_index = seg$passage_index,
        oxygen_pct = seg$oxygen_pct, recorded_duration_days = obs$passage_duration,
        selected_prediction_day = s$selection$selected_day,
        measured_initial_cells = obs$initial_cells, measured_final_cells = obs$final_cells,
        growth_scored = is.finite(obs$observed_growth), karyotype_scored = length(obs$observed_kary) > 0L,
        karyotype_cells = length(obs$observed_kary), flow_scored = present_flow,
        flow_sample_name = obs$observed_flow_sample_name,
        flow_grid_points = if (present_flow) nrow(fg) else 0L,
        counts_used_in_segment_conditioning = is.finite(obs$initial_cells) || is.finite(obs$final_cells))
      if (is.finite(obs$observed_growth)) {
        pred <- (log(s$selection$selected_live_cells) - log(s$sim$Ntot_live_obs[[1]])) / s$selection$selected_day
        sd <- fit$best_params$sigma_growth
        ll <- -log(sd) - 0.5 * log(2 * pi) - 0.5 * ((obs$observed_growth - pred) / sd)^2
        independent$growth <- c(independent$growth, ll)
        per_sample[[length(per_sample) + 1L]] <- data.frame(cohort, observation_id = pid, modality = "growth", score = ll, n_within_sample = 1L)
      }
      if (length(obs$observed_kary)) {
        prob <- smooth[as.character(obs$observed_kary)]
        prob[is.na(prob)] <- 1e-12
        ll <- mean(log(pmax(prob, 1e-12)))
        independent$karyotype <- c(independent$karyotype, ll)
        per_sample[[length(per_sample) + 1L]] <- data.frame(cohort, observation_id = pid, modality = "karyotype", score = ll, n_within_sample = length(prob))
      }
      if (present_flow) {
        grid <- as.numeric(fg$ploidy)
        density <- as.numeric(fg$density)
        keep <- is.finite(grid) & is.finite(density) & density >= 0
        stopifnot(all(keep), length(grid) > 1L)
        step <- stats::median(diff(sort(unique(grid))))
        mass <- density * step
        stopifnot(sum(mass) > 0)
        mass <- mass / sum(mass)
        pred <- vapply(grid, function(g) sum(p * stats::dnorm(g, run$grid_pre / cfg$N_UNIT, result$sigma_flow_ploidy)), numeric(1))
        stopifnot(sum(pred) > 0)
        pred <- pred / (sum(pred) * step)
        ll <- sum(mass * log(pmax(pred, 1e-12)))
        independent$flow <- c(independent$flow, ll)
        per_sample[[length(per_sample) + 1L]] <- data.frame(cohort, observation_id = pid, modality = "flow", score = ll, n_within_sample = length(grid))
      }
    }
  }
}
obs <- do.call(rbind, observations)
stopifnot(!anyDuplicated(obs$observation_id))
record("observations", "growth_count", as.integer(reference$n_growth), as.integer(sum(obs$growth_scored)), TRUE)
record("observations", "karyotype_count", as.integer(reference$n_ploidy_passages), as.integer(sum(obs$karyotype_scored)), TRUE)
record("observations", "flow_count", as.integer(reference$n_flow_samples), as.integer(sum(obs$flow_scored)), TRUE)
record("observations", "karyotyped_cells", as.integer(reference$n_kary_cells), as.integer(sum(obs$karyotype_cells)), TRUE)
ledger <- do.call(rbind, lapply(names(independent), function(n) {
  component <- switch(n, growth = "growth_loglik", karyotype = "ploidy_loglik", flow = "flow_loglik")
  value <- mean(independent[[n]])
  record("independent_score", n, result[[component]], value)
  data.frame(term = n, raw_sum = sum(independent[[n]]), denominator = length(independent[[n]]),
             mean_log_score = value, weight = 1, objective_contribution = -value, active = TRUE)
}))
for (n in c("prior", "soft_coupling", "phenotype_penalty", "failure_penalty", "death_score")) {
  ledger <- rbind(ledger, data.frame(term = n, raw_sum = 0, denominator = NA_integer_,
                  mean_log_score = 0, weight = 0, objective_contribution = 0, active = FALSE))
}
record("independent_score", "objective_sum", result$objective, sum(ledger$objective_contribution))
record("validity", "missing_growth_predictions", 0L, as.integer(result$n_growth_missing_pred), TRUE)
record("validity", "nonfinite_retained_sample_scores", 0L, as.integer(sum(!is.finite(unlist(independent)))), TRUE)
write_tsv(ledger, file.path(root, "audit", "objective_components.tsv"))
write_tsv(do.call(rbind, per_sample), file.path(out, "independent_sample_scores.tsv"))
write_tsv(obs, file.path(root, "manifests", "observation_manifest.tsv"))
write_tsv(do.call(rbind, selections), file.path(root, "audit", "selected_passage_days.tsv"))
write_tsv(do.call(rbind, boundaries), file.path(root, "audit", "conditioning_audit.tsv"))
write_tsv(do.call(rbind, distributions), file.path(out, "selected_distributions.tsv"))
tab <- read.csv(cfg$parameter_table, stringsAsFactors = FALSE)
write_tsv(tab, file.path(root, "manifests", "production_parameter_table.tsv"))
jsonlite::write_json(list(archived_cfg = fit$cfg, runtime_cfg = cfg, full_precision_parameters = fit$best_params,
                         path_relocations = list(parameter_table = list(from = original_parameter_path, to = cfg$parameter_table))),
                     file.path(root, "manifests", "runtime_model_configuration.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write_tsv(as.data.frame(installed.packages()[, c("Package", "Version", "Built")]), file.path(out, "installed_packages.tsv"))
all_checks <- do.call(rbind, checks)
write_tsv(all_checks, file.path(root, "audit", "replay_checks.tsv"))
status <- list(pass = all(all_checks$pass), checks = nrow(all_checks), failed_checks = sum(!all_checks$pass),
               optimization = FALSE, forward_evaluations = 1L, production_commit = settings$production_commit,
               passage_implementation = INVITRO_PASSAGE_IMPLEMENTATION,
               archived_objective = reference$objective, replay_objective = result$objective,
               independent_objective = sum(ledger$objective_contribution),
               forward_seconds = forward_seconds, total_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
               observations = nrow(obs), segments = length(selections), propagation_boundaries = length(boundaries),
               runtime_image = Sys.getenv("O2SD_CONTAINER_IMAGE"), slurm_job_id = Sys.getenv("SLURM_JOB_ID"),
               finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
jsonlite::write_json(status, file.path(out, "replay_status.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(status)
if (!isTRUE(status$pass)) stop("Replay checks failed; inspect audit/replay_checks.tsv")
