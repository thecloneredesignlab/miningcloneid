#!/usr/bin/env Rscript
# English diagnostic artifacts, based only on externally accepted results.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) stop("Usage: report_culture_profile_pilot.R <project> <study> [artifact-output]")
project <- normalizePath(args[1], mustWork = TRUE)
study <- normalizePath(args[2], mustWork = TRUE)
artifact_root <- if (length(args) == 3L) normalizePath(args[3], mustWork = FALSE) else study
for (folder in c("tables", "figures", "reports", "tests", "manifests")) dir.create(file.path(artifact_root, folder), recursive = TRUE, showWarnings = FALSE)
modules <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
                     c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
for (path in modules) source(path)
source(file.path(study, "scripts/culture_profile_diagnostics.R"))
settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
state_bytes <- readBin(file.path(study, "state.json"), what = "raw", n = 1e8)
state <- jsonlite::fromJSON(rawToChar(state_bytes), simplifyVector = FALSE)
state_sha256 <- digest::digest(state_bytes, algo = "sha256", serialize = FALSE)
writeBin(state_bytes, file.path(artifact_root, "manifests/report_state_snapshot.json"))
ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
dir.create(file.path(artifact_root, "figures"), showWarnings = FALSE)
library(ggplot2)
results <- list()
for (id in names(state$tasks)) if (state$tasks[[id]]$status == "accepted") {
  x <- jsonlite::read_json(file.path(study, "fits", id, "validated_result.json"), simplifyVector = TRUE)
  if (!isTRUE(x$accepted) || !isTRUE(x$launch_completed)) stop("Unaccepted endpoint in reducer")
  results[[id]] <- x
}
baseline_ids <- names(results)[vapply(names(results), function(id) state$tasks[[id]]$stage == "baseline", logical(1))]
if (!length(baseline_ids)) stop("No externally validated unrestricted baseline")
baseline_id <- baseline_ids[which.min(vapply(results[baseline_ids], function(x) x$objective, numeric(1)))]
baseline <- results[[baseline_id]]
baseline_rows <- do.call(rbind, lapply(baseline_ids, function(id) data.frame(task_id = id, role = state$tasks[[id]]$role,
  primary_family = state$tasks[[id]]$primary_family, objective = results[[id]]$objective,
  local_convergence = if (is.null(results[[id]]$optimizer_trace$local_convergence)) NA_integer_ else results[[id]]$optimizer_trace$local_convergence,
  DE_iterations = results[[id]]$optimizer_trace$iter_completed, budget_limited = results[[id]]$budget_limited)))
o2prof_write_tsv(baseline_rows, file.path(artifact_root, "tables/unrestricted_baselines.tsv"))
points <- list(); endpoint_rows <- list(); nuisance <- list(); selections <- list(); functions <- list()
growth <- list(); chromosomes <- list(); raw_kary <- list(); flow <- list(); fit_cache <- new.env(parent = emptyenv())
read_fit <- function(id) {
  if (!exists(id, envir = fit_cache, inherits = FALSE)) {
    path <- file.path(study, "fits", id, "accepted_fit.rds")
    if (!identical(o2prof_hash(path), results[[id]]$files_sha256[["accepted_fit.rds"]])) stop("Accepted fit payload changed: ", id)
    assign(id, readRDS(path), envir = fit_cache)
  }
  get(id, envir = fit_cache)
}
for (pid in names(state$points)) {
  point <- state$points[[pid]]
  ids <- unlist(point$task_ids, use.names = FALSE)
  ids <- intersect(ids, names(results))
  ids <- ids[order(vapply(results[ids], function(x) x$objective, numeric(1)))]
  if (!length(ids)) {
    points[[pid]] <- data.frame(point_id = pid, target = point$target, value = point$value,
      baseline_task_id = baseline_id, baseline_objective = baseline$objective, objective = NA_real_, delta = NA_real_,
      winner = NA_character_, accepted_starts = 0L, agreement = FALSE, active_bound_count = NA_integer_, flags = "no_verified_feasible_solution")
    next
  }
  id <- ids[1]; winner <- results[[id]]; fit <- read_fit(id)
  parameters <- read.delim(file.path(study, "fits", id, "parameter_result.tsv"), stringsAsFactors = FALSE)
  parameters$data_inactive_at_endpoint <- (parameters$parameter == "k_o_mis" & fit$best_params$p_misseg == 0) |
    (parameters$parameter == "buffer_n_exp" & fit$best_params$buffer_beta == 0) |
    (parameters$parameter %in% c("buffer_beta", "buffer_n_exp") & fit$best_params$buffer_smax == 0)
  bounds <- parameters$status == "active" & !parameters$data_inactive_at_endpoint & (parameters$at_lower | parameters$at_upper)
  agreement <- length(ids) > 1L && results[[ids[2]]]$objective - winner$objective <= settings$optimization_agreement_tolerance
  flags <- c(if (!agreement) "optimizer_agreement_unresolved", if (any(bounds)) "bound_limited",
             if (winner$budget_limited) "compute_limited", if (!isTRUE(winner$optimizer_trace$local_convergence == 0L)) "local_optimizer_unresolved")
  points[[pid]] <- data.frame(point_id = pid, target = point$target, value = point$value,
    baseline_task_id = baseline_id, baseline_objective = baseline$objective, objective = winner$objective,
    delta = winner$objective - baseline$objective, winner = id, accepted_starts = length(ids),
    agreement = agreement, active_bound_count = sum(bounds), flags = paste(flags, collapse = ";"))
  for (attempt in ids) endpoint_rows[[length(endpoint_rows) + 1L]] <- data.frame(point_id = pid, target = point$target,
    value = point$value, task_id = attempt, objective = results[[attempt]]$objective,
    delta = results[[attempt]]$objective - baseline$objective,
    growth_contribution = -results[[attempt]]$component_means$growth,
    karyotype_contribution = -results[[attempt]]$component_means$karyotype,
    flow_contribution = -results[[attempt]]$component_means$flow)
  parameters$point_id <- pid; parameters$target <- point$target; parameters$target_value <- point$value
  parameters$winning_task_id <- id
  parameters$estimated_value <- ifelse(parameters$status == "active" & !parameters$data_inactive_at_endpoint, parameters$estimate, NA_real_)
  parameters$scaled_coordinate <- NA_real_
  active <- parameters$status == "active" & !parameters$data_inactive_at_endpoint
  parameters$scaled_coordinate[active] <- (o2prof_transform(parameters$estimate[active], parameters[active, ]) -
    o2prof_transform(parameters$lower[active], parameters[active, ])) /
    (o2prof_transform(parameters$upper[active], parameters[active, ]) - o2prof_transform(parameters$lower[active], parameters[active, ]))
  nuisance[[pid]] <- parameters
  for (cohort in c("2N", "4N")) {
    run <- fit$best_components[[paste0("run_", cohort)]]
    for (segment in run$segment_results) {
      selections[[length(selections) + 1L]] <- data.frame(point_id = pid, winning_task_id = id, cohort = cohort,
        segment_id = segment$segment$segment_id, nominal_duration_days = segment$segment$duration_days,
        selected_observation_day = segment$selection$selected_day,
        parent_propagation_day = if (!is.null(segment$parent_boundary_selection)) segment$parent_boundary_selection$selected_day else NA_real_)
      probability <- segment$selection$selected_frac
      for (passage_id in segment$segment$data_ids) {
        observed <- ctx$backend$ivt_observed_passage_summary(ctx$objects$fit_data[[passage_id]])
        if (length(observed$observed_kary)) {
          display_grid <- sort(unique(c(run$grid_pre, observed$observed_kary)))
          frequency <- table(factor(observed$observed_kary, levels = display_grid)) / length(observed$observed_kary)
          K <- outer(run$grid_pre, run$grid_pre, function(a, b) dnorm(a, b, fit$best_components$sigma_kary))
          denominator <- colSums(K); denominator[denominator <= 0 | !is.finite(denominator)] <- 1
          smoothed <- as.numeric(sweep(K, 2, denominator, "/") %*% probability)
          predicted_observation <- smoothed[match(display_grid, run$grid_pre)]
          predicted_latent <- probability[match(display_grid, run$grid_pre)]
          predicted_observation[is.na(predicted_observation)] <- 0
          predicted_latent[is.na(predicted_latent)] <- 0
          chromosomes[[length(chromosomes) + 1L]] <- data.frame(point_id = pid, target = point$target,
            winning_task_id = id, cohort = cohort, passage_id = passage_id, chromosome_count = display_grid,
            observed_fraction = as.numeric(frequency), predicted_latent_fraction = predicted_latent,
            predicted_observation_fraction = predicted_observation, scoring_probability_with_floor = pmax(predicted_observation, 1e-12))
          raw_kary[[length(raw_kary) + 1L]] <- data.frame(point_id = pid, winning_task_id = id, passage_id = passage_id,
            cell_index = seq_along(observed$observed_kary), chromosome_count = observed$observed_kary)
        }
        observed_flow <- observed$observed_flow
        if (!is.null(observed_flow) && nrow(observed_flow)) {
          grid <- as.numeric(observed_flow$ploidy); step <- median(diff(sort(unique(grid))))
          predicted <- vapply(grid, function(x) sum(probability * dnorm(x, run$grid_pre / ctx$cfg$N_UNIT, fit$best_components$sigma_flow_ploidy)), numeric(1))
          predicted <- predicted / (sum(predicted) * step)
          flow[[length(flow) + 1L]] <- data.frame(point_id = pid, target = point$target, winning_task_id = id,
            cohort = cohort, passage_id = passage_id, ploidy = grid,
            observed_density = observed_flow$density / (sum(observed_flow$density) * step), predicted_density = predicted)
        }
      }
    }
  }
  g <- fit$best_components$growth_df
  g$point_id <- pid; g$target <- point$target; g$winning_task_id <- id
  growth[[pid]] <- g
  for (oxygen in c(0, 0.5, 2, 5, 20.5)) for (N in c(44L, 88L)) for (m in 0:5) {
    p <- fit$best_params
    s <- p$buffer_smax * exp(-p$buffer_beta * (2 * ctx$cfg$N_UNIT / N)^p$buffer_n_exp)
    mu <- ctx$backend$.mu_eff_of_O2(oxygen, p, N)
    saturation <- mu / (mu + p$k_o_mis)
    probability <- ctx$backend$.pmisseg_of_O2(oxygen, p, N)
    offspring <- pilotdiag_offspring_identity(N, probability, s)
    functions[[length(functions) + 1L]] <- data.frame(point_id = pid, target = point$target, target_value = point$value,
      winning_task_id = id, N = N, oxygen_pct = oxygen, error_size = m, survival = if (m == 0) 1 else s^m,
      per_error_survival = s, ploidy_survival_penalty = p$buffer_beta * (2 * ctx$cfg$N_UNIT / N)^p$buffer_n_exp,
      death_hazard_saturation = if (p$p_misseg > 0) saturation else NA_real_, induced_error_increment_before_clip = p$p_misseg * saturation,
      basal_error_probability = p$p_mis_base, death_scale_to_saturation_scale = if (p$p_misseg > 0) p$mu_hp / p$k_o_mis else NA_real_,
      p_per_chromosome = probability,
      division_rate_per_day = ctx$backend$.lambda_eff_of_O2(oxygen, p, N, O2_growth = ctx$cfg$O2_growth),
      death_rate_per_day = mu, offspring)
  }
}
combine <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
point_table <- combine(points)
for (entry in list(list(point_table, "profile_points_validated.tsv"), list(combine(endpoint_rows), "profile_endpoints.tsv"),
                  list(combine(nuisance), "nuisance_paths.tsv"), list(combine(selections), "selected_passage_days.tsv"),
                  list(combine(functions), "state_function_values.tsv"), list(combine(growth), "growth_predictions.tsv"),
                  list(combine(chromosomes), "chromosome_distributions.tsv"), list(combine(raw_kary), "observed_karyotype_cells.tsv"),
                  list(combine(flow), "flow_distributions.tsv")))
  o2prof_write_tsv(entry[[1]], file.path(artifact_root, "tables", entry[[2]]))
if (nrow(point_table)) {
  for (target in c("p_misseg", "buffer_beta")) {
    p <- point_table[point_table$target == target, ]; p <- p[order(p$value), ]
    display_scale <- if (target == "p_misseg") 1e-4 else 1e-3
    p$coordinate <- log1p(p$value / display_scale)
    p$line_group <- cumsum(!p$agreement | !is.finite(p$delta))
    endpoints <- combine(endpoint_rows); endpoints <- endpoints[endpoints$target == target, ]
    endpoints$coordinate <- log1p(endpoints$value / display_scale)
    tick_rows <- 1L
    for (i in seq_len(nrow(p))[-1]) if (p$coordinate[i] - p$coordinate[tail(tick_rows, 1)] >= .35) tick_rows <- c(tick_rows, i)
    if (!nrow(p) %in% tick_rows) {
      tick_rows <- tick_rows[p$coordinate[nrow(p)] - p$coordinate[tick_rows] >= .35]
      tick_rows <- c(tick_rows, nrow(p))
    }
    plot_prefix <- if (state$status == "pilot_fitting_finished") "Culture" else if (state$status %in% c("baseline_optimizer_unresolved", "baseline_repair_budget_exhausted")) "Incomplete culture" else "Interim culture"
    chart <- ggplot() + geom_point(data = endpoints, aes(coordinate, delta), colour = "grey70", size = 1, alpha = .5) +
      geom_line(data = p[p$agreement & is.finite(p$delta), ], aes(coordinate, delta, group = line_group), colour = "#194B7A") +
      geom_point(data = p[is.finite(p$delta), ], aes(coordinate, delta, shape = agreement), colour = "#194B7A", size = 2.5) +
      geom_point(data = data.frame(coordinate = log1p(read_fit(baseline_id)$best_params[[target]] / display_scale)),
                 aes(coordinate, y = 0), shape = 8, colour = "#B35B16", size = 3) +
      scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "Start agreement") +
      scale_x_continuous(breaks = p$coordinate[tick_rows], labels = formatC(p$value[tick_rows], digits = 3, format = "g")) +
      labs(x = paste0(target, " (zero-capable spacing)"), y = "Best verified objective minus unrestricted baseline",
           title = paste(plot_prefix, "profile objective:", target),
           subtitle = "Grey: endpoints; star: unrestricted candidate at its parameter value; no confidence thresholds") + theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(artifact_root, "figures", paste0(target, "_profile.png")), chart, width = 8, height = 5, dpi = 220)
    if (any(is.finite(p$delta))) {
      zoom <- chart + coord_cartesian(ylim = range(c(0, p$delta[is.finite(p$delta)]))) +
        labs(subtitle = "Display zoom to verified minima; all endpoints retained in the overview; no confidence thresholds")
      ggsave(file.path(artifact_root, "figures", paste0(target, "_profile_minima.png")), zoom, width = 8, height = 5, dpi = 220)
    }
    paths <- combine(nuisance); paths <- paths[paths$target == target & paths$parameter != target, ]
    paths$scan_value <- factor(paths$target_value, levels = p$value, labels = formatC(p$value, digits = 3, format = "g"))
    chart <- ggplot(paths, aes(scan_value, parameter, fill = scaled_coordinate)) + geom_tile() +
      scale_x_discrete(drop = FALSE) +
      scale_fill_gradient(low = "#EFF3F8", high = "#194B7A", limits = c(0, 1), na.value = "grey80", name = "Fraction of\ntransformed range") +
      labs(x = paste0(target, " (fixed scan values)"), y = NULL, title = "Nuisance compensation along verified minima",
           subtitle = "Grey: data-inactive; blank columns: no verified result; bounds remain part of the study") + theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(artifact_root, "figures", paste0(target, "_compensation.png")), chart, width = 8, height = 7, dpi = 220)
  }
}
# Dense flux is exported for the baseline plus null, high-end and largest-rise representatives.
representatives <- baseline_id
if (nrow(point_table)) for (target in c("p_misseg", "buffer_beta")) {
  p <- point_table[point_table$target == target & is.finite(point_table$objective), ]
  if (nrow(p)) representatives <- c(representatives, p$winner[which.min(p$value)], p$winner[which.max(p$value)], p$winner[which.max(p$delta)])
}
representatives <- unique(representatives)
representative_points <- if (nrow(point_table)) point_table$point_id[point_table$winner %in% representatives] else character()
responses <- list(); survivals <- list()
oxygen_grid <- sort(unique(c(seq(0, 5, by = .025), exp(seq(log(.001), log(20.5), length.out = 100)), 20.5)))
for (id in representatives) {
  fit <- read_fit(id); p <- fit$best_params
  selected_point <- if (nrow(point_table)) point_table[!is.na(point_table$winner) & point_table$winner == id, , drop = FALSE] else data.frame()
  if (id != baseline_id && !nrow(selected_point)) stop("Representative has no verified profile point: ", id)
  label <- if (id == baseline_id) paste0("Unrestricted (", id, ")") else
    paste0(selected_point$target[1], " = ", formatC(selected_point$value[1], digits = 3, format = "g"), " (", id, ")")
  for (N in c(44L, 88L)) {
    mu <- ctx$backend$.mu_eff_of_O2(oxygen_grid, p, N)
    probability <- ctx$backend$.pmisseg_of_O2(oxygen_grid, p, N)
    s <- p$buffer_smax * exp(-p$buffer_beta * (2 * ctx$cfg$N_UNIT / N)^p$buffer_n_exp)
    responses[[length(responses) + 1L]] <- data.frame(task_id = id, fit_label = label, N = N, oxygen_pct = oxygen_grid,
      p_per_chromosome = probability, death_rate_per_day = mu,
      division_rate_per_day = ctx$backend$.lambda_eff_of_O2(oxygen_grid, p, N, O2_growth = ctx$cfg$O2_growth),
      induced_error_increment_before_clip = p$p_misseg * mu / (mu + p$k_o_mis),
      pilotdiag_offspring_identity(N, probability, s))
  }
  N <- ctx$cfg$N_MIN:ctx$cfg$N_MAX
  s <- p$buffer_smax * exp(-p$buffer_beta * (2 * ctx$cfg$N_UNIT / N)^p$buffer_n_exp)
  for (m in 0:5) survivals[[length(survivals) + 1L]] <- data.frame(task_id = id, fit_label = label,
    N = N, error_size = m, per_error_survival = s, survival = if (m == 0) rep(1, length(N)) else s^m)
}
response_table <- combine(responses); survival_table <- combine(survivals)
o2prof_write_tsv(response_table, file.path(artifact_root, "tables/effective_response_functions.tsv"))
o2prof_write_tsv(survival_table, file.path(artifact_root, "tables/survival_functions.tsv"))
response_long <- do.call(rbind, lapply(c("p_per_chromosome", "death_rate_per_day", "division_rate_per_day"), function(quantity) {
  frame <- response_table[, c("task_id", "fit_label", "N", "oxygen_pct")]
  frame$quantity <- quantity; frame$value <- response_table[[quantity]]; frame
}))
response_long$quantity <- factor(response_long$quantity,
  levels = c("p_per_chromosome", "death_rate_per_day", "division_rate_per_day"),
  labels = c("Error probability\n(per chromosome)", "Death hazard\n(per day)", "Division rate\n(per day)"))
chart <- ggplot(response_long, aes(log1p(oxygen_pct / .025), value, colour = fit_label)) + geom_line(linewidth = .6) +
  facet_grid(quantity ~ N, scales = "free_y", labeller = labeller(quantity = label_value, N = label_both)) +
  scale_x_continuous(breaks = log1p(c(0, .1, .5, 2, 5, 20.5) / .025), labels = c(0, .1, .5, 2, 5, 20.5)) +
  theme_bw(base_size = 9) + labs(x = "Assigned oxygen (%) with zero-capable spacing", y = "Value (units shown in row labels)",
    colour = "Verified candidate", title = "Effective chromosome-error, death and division functions",
    subtitle = "Matched parental chromosome counts; candidate traces, not confidence bands")
ggsave(file.path(artifact_root, "figures/effective_response_functions.png"), chart, width = 13, height = 9, dpi = 220)
offspring_long <- do.call(rbind, lapply(c("effective_error_probability_conditional_on_viability", "biological_viable_daughters_per_ordinary_division"), function(quantity) {
  frame <- response_table[, c("task_id", "fit_label", "N", "oxygen_pct")]
  frame$quantity <- quantity; frame$value <- response_table[[quantity]]; frame
}))
offspring_long$quantity <- factor(offspring_long$quantity,
  levels = c("effective_error_probability_conditional_on_viability", "biological_viable_daughters_per_ordinary_division"),
  labels = c("Conditional error probability\n(per chromosome)", "Biologically viable daughters\n(per ordinary division)"))
chart <- ggplot(offspring_long, aes(log1p(oxygen_pct / .025), value, colour = fit_label)) + geom_line(linewidth = .6) +
  facet_grid(quantity ~ N, scales = "free_y", labeller = labeller(quantity = label_value, N = label_both)) +
  scale_x_continuous(breaks = log1p(c(0, .1, .5, 2, 5, 20.5) / .025), labels = c(0, .1, .5, 2, 5, 20.5)) +
  theme_bw(base_size = 9) + labs(x = "Assigned oxygen (%) with zero-capable spacing", y = "Value (units shown in row labels)",
    colour = "Verified candidate", title = "Effective ordinary-daughter functions",
    subtitle = "Before numerical truncation/grid loss; diagnostic traces, not confidence bands")
ggsave(file.path(artifact_root, "figures/effective_daughter_functions.png"), chart, width = 13, height = 7, dpi = 220)
chart <- ggplot(survival_table, aes(N, survival, colour = fit_label)) + geom_line(linewidth = .6) +
  facet_wrap(~ error_size, ncol = 2, labeller = label_both) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_vline(xintercept = c(44, 88), linetype = "dotted", colour = "grey70") + theme_bw(base_size = 9) +
  labs(x = "Parental chromosome count", y = "Survival conditional on error size", colour = "Verified candidate",
       title = "Survival functions across the production chromosome grid",
       subtitle = "Error size zero has survival one; candidate traces, not confidence bands")
ggsave(file.path(artifact_root, "figures/survival_functions.png"), chart, width = 13, height = 9, dpi = 220)
g <- combine(growth)
if (nrow(g)) {
  g <- g[g$point_id %in% representative_points, ]
  chart <- ggplot(g, aes(observed_growth, predicted_growth_rate, colour = factor(oxygen_pct), shape = cohort)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70") + geom_point(size = 1.2, alpha = .7) +
    facet_wrap(~ point_id, scales = "free") + theme_bw(base_size = 9) +
    labs(x = "Observed interval growth rate (per day)", y = "Predicted interval growth rate (per day)",
         colour = "Assigned oxygen (%)", title = "Growth predictions at representative restricted fits")
  ggsave(file.path(artifact_root, "figures/growth_predictions.png"), chart, width = 10, height = 8, dpi = 220)
}
for (pid in representative_points) {
  k <- combine(chromosomes); k <- k[k$point_id == pid, ]
  if (nrow(k)) {
    chart <- ggplot(k, aes(chromosome_count)) + geom_col(aes(y = observed_fraction), fill = "grey75", width = .9) +
      geom_line(aes(y = predicted_observation_fraction), colour = "#194B7A", linewidth = .5) +
      facet_wrap(~ passage_id, ncol = 3, scales = "free_y", labeller = label_wrap_gen(28)) + theme_bw(base_size = 8) +
      labs(x = "Chromosome count", y = "Sample fraction", title = paste("Full chromosome distributions:", pid),
           subtitle = "Grey bars: observations; blue: predicted observation distribution")
    ggsave(file.path(artifact_root, "figures", paste0(pid, "_chromosomes.png")), chart, width = 10, height = 10, dpi = 200)
  }
  f <- combine(flow); f <- f[f$point_id == pid, ]
  if (nrow(f)) {
    chart <- ggplot(f, aes(ploidy)) + geom_line(aes(y = observed_density), colour = "grey55") +
      geom_line(aes(y = predicted_density), colour = "#194B7A") +
      facet_wrap(~ passage_id, ncol = 4, scales = "free_y", labeller = label_wrap_gen(25)) + theme_bw(base_size = 8) +
      labs(x = "Ploidy units (chromosome count / 22)", y = "Normalized density", title = paste("Full flow distributions:", pid),
           subtitle = "Grey: observations; blue: predicted density")
    ggsave(file.path(artifact_root, "figures", paste0(pid, "_flow.png")), chart, width = 12, height = 10, dpi = 200)
  }
}
flux_summaries <- list(); flux_series <- list(); kernel_tables <- list(); dense_checks <- list()
for (id in representatives) {
  fit <- read_fit(id); cfg <- fit$cfg
  kernels <- new.env(parent = emptyenv())
  for (cohort in c("2N", "4N")) for (segment in fit$best_components[[paste0("run_", cohort)]]$segment_results) {
    key <- as.character(segment$segment$oxygen_pct)
    if (!exists(key, kernels, inherits = FALSE)) {
      kernel <- pilotdiag_rates(ctx, fit$best_params, cfg, segment$segment$oxygen_pct)
      assign(key, kernel, kernels)
      table <- kernel$rates; table$task_id <- id; kernel_tables[[length(kernel_tables) + 1L]] <- table
    }
    audited <- pilotdiag_dense_flux(ctx, fit, cohort, segment, get(key, kernels))
    audited$flux$task_id <- id; audited$summary$task_id <- id
    flux_series[[length(flux_series) + 1L]] <- audited$flux
    flux_summaries[[length(flux_summaries) + 1L]] <- audited$summary
    dense_checks[[length(dense_checks) + 1L]] <- data.frame(task_id = id, cohort = cohort, segment_id = segment$segment$segment_id,
      maximum_saved_time_state_difference = audited$daily_replay_maximum_difference)
  }
}
o2prof_write_tsv(combine(flux_summaries), file.path(artifact_root, "tables/passage_flux_summaries.tsv"))
o2prof_write_tsv(combine(flux_series), file.path(artifact_root, "tables/daughter_fluxes.tsv"))
o2prof_write_tsv(combine(kernel_tables), file.path(artifact_root, "tables/state_daughter_rates.tsv"))
o2prof_write_tsv(combine(dense_checks), file.path(artifact_root, "tests/final_dense_replay_checks.tsv"))
o2pilot_json(list(pass = TRUE, baseline_task_id = baseline_id, baseline_objective = baseline$objective,
                  execution_status = state$status, report_state_sha256 = state_sha256,
                  baseline_search_gate_passed = isTRUE(state$baseline_gate$pass_gate),
                  point_count = nrow(point_table), points_without_verified_solution = if (nrow(point_table)) sum(!is.finite(point_table$objective)) else 0,
                  points_without_start_agreement = if (nrow(point_table)) sum(!point_table$agreement) else 0,
                  flux_representatives = representatives, flux_segment_count = length(flux_summaries),
                  inference_calibration = "absent", output_language = "English"), file.path(artifact_root, "reports/diagnostic_artifacts_status.json"))
cat("English profile, compensation, prediction and event-accounting artifacts written.\n")
