# Static reporting from saved, independently validated artifacts. No simulation.
report_dir <- file.path(study, "reports"); dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
cs <- read_table(file.path(study, "tables/recovery_candidates.tsv"))
for (i in seq_len(nrow(cs))) if (cs$mode[i] != "reference") {
  d <- dirname(file.path(project, cs$source_fit[i])); payload <- read_table(file.path(d, "payload_manifest.tsv"))
  if (!all(vapply(file.path(d, payload$path), o2prof_hash, character(1)) == payload$sha256)) stop("Final payload fingerprint mismatch")
}
gate <- read_kv(file.path(study, "tables/recovery_gate.tsv"))
settings <- read_kv(file.path(study, "settings.tsv"))
quality <- list()
for (i in seq_len(nrow(cs))) if (cs$mode[i] != "reference") {
  d <- dirname(file.path(project, cs$source_fit[i])); local <- file.path(d, "local_search.rds")
  trace <- list.files(file.path(d, "best_candidates/improvements"), pattern = "\\.rds$")
  if (file.exists(local)) {
    l <- readRDS(local); sb <- readRDS(file.path(d, "sbplx.rds")); method <- read_kv(file.path(d, "method_sequence.tsv"))
    quality[[length(quality) + 1L]] <- data.frame(task_id = cs$task_id[i], origin = cs$origin[i],
      family = cs$family[i], sbplx_status = sb$status, sbplx_message = sb$message,
      full_objective_calls = l$evaluations, local_stable = l$tested_directional_stability,
      overall_budget_limited = l$budget_limited, improvement_records = length(trace),
      improvement_archive = if (length(trace)) "append_only_full_parameter_component_history" else "historical_diagnostic_endpoint_only", stringsAsFactors = FALSE)
  }
}
if (length(quality)) o2prof_write_tsv(do.call(rbind, quality), file.path(report_dir, "local_optimization_quality.tsv"))
cs <- cs[order(cs$objective, cs$task_id), ]
cs <- cs[!duplicated(cs$source_fit), ]
selected <- cs[!duplicated(cs$origin), ]; selected <- head(selected, 4L)
all_scores <- list(); all_params <- list(); days <- list(); growth <- list(); distributions <- list(); flows <- list()
for (i in seq_len(nrow(selected))) {
  z <- selected[i, ]; fit <- readRDS(file.path(project, z$source_fit)); comp <- fit$best_components
  all_scores[[i]] <- data.frame(task_id = z$task_id, origin = z$origin, t(score_components(comp)), check.names = FALSE)
  map <- original_map(ctx); tab <- map$table
  tab$estimate <- vapply(tab$parameter, function(n) fit$best_params[[n]], numeric(1)); tab$task_id <- z$task_id; tab$origin <- z$origin
  tab$normalized_coordinate <- o2review_unit_spec(map, fit$best_params)$unit; all_params[[i]] <- tab
  g <- comp$growth_df; g$task_id <- z$task_id; growth[[i]] <- g
  f <- comp$flow_overlay_df; f$task_id <- z$task_id; flows[[i]] <- f
  for (cohort in c("2N", "4N")) {
    run <- comp[[paste0("run_", cohort)]]
    for (j in seq_along(run$segment_results)) {
      seg <- run$segment_results[[j]]
      days[[length(days) + 1L]] <- data.frame(task_id = z$task_id, origin = z$origin, cohort = cohort,
        segment_order = j, segment_id = seg$segment$segment_id, oxygen_pct = seg$segment$oxygen_pct,
        selected_day = seg$selection$selected_day, selected_live_cells = seg$selection$selected_live_cells)
      distributions[[length(distributions) + 1L]] <- data.frame(task_id = z$task_id, origin = z$origin, cohort = cohort,
        segment_order = j, segment_id = seg$segment$segment_id, chromosome_number = run$grid_pre,
        live_fraction = seg$selection$selected_frac, propagated_cells = seg$selection$reseeded_state)
    }
  }
}
score <- do.call(rbind, all_scores); parameters <- do.call(rbind, all_params); day <- do.call(rbind, days)
o2prof_write_tsv(score, file.path(report_dir, "representative_components.tsv"))
o2prof_write_tsv(parameters, file.path(report_dir, "representative_parameters.tsv"))
o2prof_write_tsv(day, file.path(report_dir, "representative_selected_days.tsv"))
o2prof_write_tsv(do.call(rbind, growth), file.path(report_dir, "representative_growth.tsv"))
o2prof_write_tsv(do.call(rbind, distributions), file.path(report_dir, "representative_full_distributions.tsv"))
o2prof_write_tsv(do.call(rbind, flows), file.path(report_dir, "representative_flow_overlays.tsv"))
# Assemble historical global endpoints and current global/complete-chain scores.
ledger <- read_table(file.path(study, "manifests/historical_origins.tsv")); chains <- list()
origins <- unique(cs$origin)
for (origin in origins) {
  r <- cs[cs$origin == origin, ]; local <- r[r$mode == "local", ]; global <- r[r$mode == "global", ]
  old <- ledger[ledger$task_id == origin & ledger$pure_independent %in% c(TRUE, "True", "TRUE"), ]
  raw <- if (nrow(global)) min(global$objective) else if (nrow(old)) old$objective[1] else NA_real_
  q <- if (nrow(local)) min(local$objective) else NA_real_
  chains[[length(chains) + 1L]] <- data.frame(origin = origin, family = r$family[1], global_objective = raw,
    complete_chain_objective = q, global_gap = raw - gate$q_star, chain_gap = q - gate$q_star,
    prospective = r$family[1] %in% c("A", "B"), local_stable = if (nrow(local)) local$local_stable[which.min(local$objective)] else FALSE)
}
chain <- do.call(rbind, chains); o2prof_write_tsv(chain, file.path(report_dir, "global_vs_complete_chain.tsv"))
pdf(file.path(report_dir, "baseline_diagnostics.pdf"), width = 11, height = 8.5, useDingbats = FALSE)
par(mar = c(6, 5, 4, 1))
x <- seq_len(nrow(chain)); y <- range(c(chain$global_objective, chain$complete_chain_objective, gate$q_star), finite = TRUE)
plot(x, chain$global_objective, pch = 1, xaxt = "n", ylim = y, xlab = "", ylab = "Data-only objective Q", main = "Global endpoints versus complete optimization chains")
points(x, chain$complete_chain_objective, pch = 19, col = "steelblue")
segments(x, chain$global_objective, x, chain$complete_chain_objective, col = "grey60")
axis(1, at = x, labels = chain$origin, las = 2); abline(h = gate$q_star + .001, col = "firebrick", lty = 2)
legend("topright", c("Global endpoint", "Complete chain", "Final reference + 0.001"), pch = c(1, 19, NA), lty = c(NA, NA, 2), col = c("black", "steelblue", "firebrick"), bty = "n")
par(mar = c(9, 4, 4, 1)); colors <- seq_len(nrow(selected))
mat <- do.call(cbind, lapply(all_params, function(z) z$normalized_coordinate))
matplot(seq_len(nrow(mat)), mat, type = "b", pch = seq_len(ncol(mat)), col = colors, lty = 1, xaxt = "n", ylim = c(0, 1), xlab = "", ylab = "Fraction of transformed parameter range", main = "Representative parameter regions; no posterior interpretation")
axis(1, at = seq_len(nrow(mat)), labels = all_params[[1]]$parameter, las = 2, cex.axis = .7)
legend("topright", selected$origin, col = colors, pch = seq_len(ncol(mat)), lty = 1, bty = "n")
par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
for (cohort in c("2N", "4N")) {
  curves <- lapply(seq_len(nrow(selected)), function(i) {
    fit <- readRDS(file.path(project, selected$source_fit[i])); run <- fit$best_components[[paste0("run_", cohort)]]
    list(grid = run$grid_pre, fraction = tail(run$segment_results, 1)[[1]]$selection$selected_frac)
  })
  matplot(curves[[1]]$grid, do.call(cbind, lapply(curves, function(z) z$fraction)), type = "l", col = colors, lty = 1,
    xlab = "Chromosome number", ylab = "Live-cell fraction", main = paste(cohort, "last saved segment"))
  legend("topright", selected$origin, col = colors, lty = 1, bty = "n", cex = .8)
}
par(mfrow = c(1, 1), mar = c(5, 5, 4, 1))
g <- do.call(rbind, growth); plot(g$observed_growth, g$predicted_growth_rate, pch = 16, cex = .6, col = match(g$task_id, selected$task_id),
  xlab = "Observed growth rate", ylab = "Saved predicted growth rate", main = "Representative saved growth predictions")
abline(0, 1, lty = 2)
dev.off()
status <- if (isTRUE(gate$pass_gate)) "PASSED" else "UNRESOLVED"
lines <- c("# Unrestricted culture baseline recovery", "", paste("Prospective complete-chain recovery:", status), "",
  sprintf("The lowest independently verified feasible objective is %.17g. Eligible prospective origins: %d. Represented families: %s. Strong final-repeat criterion: %s.", gate$q_star, gate$eligible_origins, gate$families, gate$strong_repeat), "",
  "## Workflow and numerical validation", "",
  "The executed workflow uses R, Python, shell and the pinned production-v2 C++ model. It creates no Java, JavaScript or JSON artifacts or runtime dependencies. Typed TSV transactions and RDS carry state. Complete payloads require a fresh-process score, component, observation, selected-day and propagated-state check. Numerical tolerances were not relaxed. Cache build/hit counters are compared separately as process-history telemetry; they are not predictions or propagated states.", "",
  "Production DE retains its existing graded protocol feasibility guidance. Direct local evaluation keeps Inf for infeasible probes. No new finite penalty, parameter-domain extension, prior or frozen selected-day approximation was introduced. The model, full culture data, modality weights and protocol rules remain frozen.", "",
  "## What was tested", "",
  "Six historical endpoint-specific local jobs diagnose the old global-only recovery gap. New A/B chains use predeclared independent populations, native DE global search and separately allocated native Sbplx plus bounded pattern/boundary checks. Independent-origin count, local stability, native iteration limits and fresh-process reproducibility are recorded separately. The old global-only recovery failure is preserved; retrospective polishing does not retroactively pass it.", "",
  "The first six historical diagnostics retain evaluation scores, selected-day records and verified final vectors; their intermediate improvement vectors were not archived by the initial engineering version. This limitation is preserved explicitly. Prospective confirmation requires append-only full parameter, component-score and selected-day records. See local_optimization_quality.tsv for archive coverage, native termination status and tested local stability.", "",
  "See global_vs_complete_chain.tsv for every origin, representative_parameters.tsv for parameter regions, and the full fit directories for native starts, evaluation logs, method sequence, local sweeps, boundary probes, validation receipts and source fingerprints. These are optimization diagnostics; the endpoints are not posterior samples or parameter confidence intervals.", "",
  "## Saved predictions and interpretation", "",
  "baseline_diagnostics.pdf shows endpoint scores, transformed parameter regions, saved full chromosome distributions and observed-versus-predicted growth. Full representative growth, flow overlays, distributions and selected days are exported as separate keyed tables. Similar scores with different coordinates indicate remaining compensation in this declared search and data setting; they do not prove structural nonidentifiability.", "",
  "## Resource closure and scope", "",
  "All allocated tasks use QOS xxlarge. There is no application concurrency cap. Tables/jobs.tsv preserves Slurm state, exit codes, reservations and allocated core-hours separately from internal acceptance. Final accounting must be refreshed after the controller and report allocations close. The new hard budget is 1440 core-hours, with a 120-core-hour non-fit sub-cap and declared attempt limits.", "",
  "No constrained profiles, tumor fitting, COMBOS, measurement deletion or statistical calibration was started. This baseline-only phase ends with this report; later work requires a separate scope decision.")
writeLines(lines, file.path(study, "REPORT.md"))
assert_runtime(report_dir); writeLines("Completed static report without model optimization", file.path(report_dir, "COMPLETE"))
