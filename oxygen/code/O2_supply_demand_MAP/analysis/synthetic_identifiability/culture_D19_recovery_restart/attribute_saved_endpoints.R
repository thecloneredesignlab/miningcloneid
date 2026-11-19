# Archived-state attribution. This entry point never propagates or fits the model.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
study <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
stopifnot(!file.exists(out))

read_table <- function(rel) {
  read.delim(file.path(study, rel), check.names = FALSE, stringsAsFactors = FALSE)
}
index <- read_table('local_collection_manifest.tsv')
delivery_index <- read_table('reports/post_delivery_manifest.tsv')
index <- rbind(index, delivery_index[!delivery_index$path %in% index$path, names(index)])
used <- character()
sha256 <- function(path) {
  sha <- Sys.which('sha256sum')
  if (nzchar(sha)) {
    command <- sha
    flags <- character()
  } else {
    command <- Sys.which('shasum')
    flags <- c('-a', '256')
    stopifnot(nzchar(command))
  }
  result <- system2(command, c(flags, shQuote(path)), stdout = TRUE)
  stopifnot(length(result) == 1L, is.null(attr(result, 'status')))
  hash <- sub(' .*', '', result)
  stopifnot(grepl('^[a-f0-9]{64}$', hash))
  hash
}
verified_path <- function(rel) {
  row <- index[index$path == rel, , drop = FALSE]
  path <- file.path(study, rel)
  stopifnot(nrow(row) == 1L, file.info(path)$size == row$bytes, sha256(path) == row$sha256)
  used <<- unique(c(used, rel))
  path
}
read_archive <- function(rel) readRDS(verified_path(rel))
read_verified <- function(rel) read.delim(verified_path(rel), check.names = FALSE)
lock <- read_verified('reports/selection_lock.tsv')
gate <- read_verified('reports/final_gate.tsv')
stopifnot(nrow(lock) == 1L, nrow(gate) == 1L, lock$status == 'DET_OPTIMIZER_UNRESOLVED',
          gate$status == lock$status, gate$active_owned_jobs == 0L,
          identical(trimws(readLines(verified_path('controller/FINISHED'))), lock$status))
best <- read_verified('reports/best_found_gate_metrics.tsv')
quality <- read_verified('reports/verified_candidates.tsv')
stopifnot(identical(sort(best$origin), paste0('D_', LETTERS[1:6])),
          all(best$fresh_replay_verified), all(!best$prediction_pass),
          identical(best$key[which.min(best$D)], lock$best_found_key))

case <- read_archive('inputs/public/D19.rds')
contract <- read_archive('inputs/public/contract.rds')
truth <- read_archive('inputs/private/truth.rds')
reference <- read_archive('inputs/private/reference.rds')$pred
growth_scale <- case$known_sigma_growth
grid <- 22:154
stopifnot(identical(growth_scale, truth$sigma_growth), growth_scale > 0,
          identical(case$data$growth, reference$growth),
          identical(case$data$prob, reference$kary),
          length(case$data$growth) == 114L, length(case$data$prob) == 12L)

# Evaluate only the published observer function definition from the frozen source.
observer_source <- verified_path('inputs/production/o2_supply_demand_map_invitro_objective_utils.R')
observer_definition <- Filter(function(expr) {
  is.call(expr) && identical(expr[[1L]], as.name('<-')) &&
    identical(expr[[2L]], as.name('ivt_smooth_kary_distribution'))
}, as.list(parse(observer_source)))
stopifnot(length(observer_definition) == 1L)
observer_environment <- new.env(parent = baseenv())
eval(observer_definition[[1L]], observer_environment)
observe <- function(p, width) observer_environment$ivt_smooth_kary_distribution(grid, p, width)

segments <- do.call(c, unname(lapply(contract$design$adapters, function(x) x$segments)))
segment_key <- vapply(segments, function(x) paste(x$cohort, x$segment_id), character(1))
edge_key <- with(contract$edges, paste(cohort, child))
stopifnot(length(segments) == 122L, !anyDuplicated(segment_key),
          setequal(edge_key, segment_key))
edges <- contract$edges[match(segment_key, edge_key), , drop = FALSE]
parent_index <- match(with(edges, paste(cohort, parent)), segment_key)
stopifnot(sum(!is.na(parent_index)) == 120L,
          all(is.na(parent_index) == (edges$kind == 'root')))
passage_ids <- lapply(segments, function(x) x$data_ids)
passage_lookup <- setNames(rep(seq_along(segments), lengths(passage_ids)),
                           unlist(passage_ids, use.names = FALSE))
stopifnot(!anyDuplicated(names(passage_lookup)))

normalized <- function(x) {
  stopifnot(length(x) == length(grid), all(is.finite(x)), all(x >= 0), sum(x) > 0)
  as.numeric(x / sum(x))
}
same <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
tv <- function(a, b) sum(abs(a - b)) / 2
hellinger2 <- function(a, b) sum((sqrt(a) - sqrt(b))^2) / 2

validate_snapshot <- function(pred) {
  selection <- pred$selection
  key <- with(selection, paste(cohort, segment_id))
  stopifnot(!anyDuplicated(key), setequal(key, segment_key),
            setequal(names(pred$states), segment_key))
  selection <- selection[match(segment_key, key), , drop = FALSE]
  snapshots <- pred$states[segment_key]
  totals <- lapply(snapshots, rowSums)
  reseed_mode <- rep('root', length(segments))
  choice_margin_cells <- numeric(length(segments))
  reseed_state_error <- 0
  for (i in seq_along(segments)) {
    days <- segments[[i]]$obs_days_local
    n <- totals[[i]]
    stopifnot(ncol(snapshots[[i]]) == length(grid), length(days) == length(n),
              abs(n[1] / edges$initial_cells_effective[i] - 1) <= 1e-10,
              selection$day[i] == days[selection$index[i]])
    candidates <- which(days > 0)
    if (!length(candidates)) candidates <- seq_along(days)
    target <- edges$target_effective[i]
    if (is.finite(target) && target > 0) {
      gaps <- abs(n[candidates] - target)
      ordered <- order(gaps, candidates)
      expected <- candidates[ordered[1L]]
      choice_margin_cells[i] <- if (length(ordered) > 1L) diff(gaps[ordered[1:2]]) else NA_real_
    } else {
      expected <- tail(candidates, 1L)
      choice_margin_cells[i] <- NA_real_
    }
    stopifnot(selection$index[i] == expected)
    p <- parent_index[i]
    if (is.na(p)) next
    parent_days <- segments[[p]]$obs_days_local
    parent_total <- totals[[p]]
    selected <- selection$index[p]
    required <- edges$initial_cells_effective[i]
    if (parent_total[selected] >= required) {
      expected_parent <- selected
      reseed_mode[i] <- 'target_selected_parent_day'
    } else {
      eligible <- which(parent_days > 0)
      if (!length(eligible)) eligible <- seq_along(parent_days)
      eligible <- eligible[parent_total[eligible] >= required]
      stopifnot(length(eligible) > 0L)
      expected_parent <- eligible[order(parent_total[eligible] - required, eligible)[1L]]
      reseed_mode[i] <- 'smallest_eligible_surplus'
    }
    stopifnot(selection$parent_index[i] == expected_parent,
              selection$parent_day[i] == parent_days[expected_parent])
    actual <- normalized(snapshots[[i]][1, ])
    expected_state <- normalized(snapshots[[p]][expected_parent, ])
    reseed_state_error <- max(reseed_state_error, tv(actual, expected_state))
  }
  stopifnot(reseed_state_error <= 1e-10)
  list(selection = selection, snapshots = snapshots, totals = totals,
       reseed_mode = reseed_mode, choice_margin_cells = choice_margin_cells,
       reseed_state_error = reseed_state_error)
}

ancestry <- function(i) {
  result <- integer()
  while (!is.na(i)) {
    stopifnot(!i %in% result)
    result <- c(i, result)
    i <- parent_index[i]
  }
  result
}
growth_at <- function(snapshot, segment, row) {
  total <- snapshot$totals[[segment]]
  day <- segments[[segment]]$obs_days_local[row]
  stopifnot(day > 0, total[1] > 0, total[row] > 0)
  log(total[row] / total[1]) / day
}

base <- validate_snapshot(reference)
per_segment <- list()
per_growth <- list()
per_kary <- list()
summaries <- list()
validation <- list()
search_stages <- list()

for (j in seq_len(nrow(best))) {
  b <- best[j, ]
  row <- quality[quality$origin == b$origin & quality$key == b$key, , drop = FALSE]
  stopifnot(nrow(row) == 1L, row$verified, row$D == b$D)
  saved <- read_archive(row$state_file)
  candidate <- saved$candidate
  pred <- saved$pred
  stopifnot(identical(candidate$key, b$key), identical(candidate$origin, b$origin),
            identical(candidate$params$sigma_growth, growth_scale),
            identical(names(pred$growth), names(case$data$growth)),
            identical(names(pred$kary), names(case$data$prob)))
  stages <- c('bench', 'de6', 'local', 'pattern', 'replay2b')
  search_stages[[j]] <- do.call(rbind, lapply(stages, function(stage) {
    phase <- paste0(b$origin, '_', stage)
    state <- read_archive(file.path('verification', phase, 'state.rds'))
    x <- state$best
    stopifnot(identical(state$origin, b$origin), identical(state$stage, stage),
              is.finite(x$score), identical(x$params$sigma_growth, growth_scale))
    data.frame(origin = b$origin, stage = stage, key = x$key, D = x$score,
               D_G = x$loss$growth, D_K = x$loss$kary, sigma_kary = x$params$sigma_kary,
               init_sd_2N = x$params$init_sd_2N, init_sd_4N = x$params$init_sd_4N)
  }))
  stopifnot(all(diff(search_stages[[j]]$D) <= 0),
            identical(tail(search_stages[[j]]$key, 1L), b$key))
  fitted <- validate_snapshot(pred)
  obs_shift <- !same(fitted$selection$index, base$selection$index) |
    !same(fitted$selection$day, base$selection$day)
  reset_shift <- !same(fitted$selection$parent_index, base$selection$parent_index) |
    !same(fitted$selection$parent_day, base$selection$parent_day)
  mode_shift <- fitted$reseed_mode != base$reseed_mode
  path_reset <- vapply(seq_along(segments), function(i) sum(reset_shift[ancestry(i)]), integer(1))
  path_mode <- vapply(seq_along(segments), function(i) sum(mode_shift[ancestry(i)]), integer(1))
  first_reset <- vapply(seq_along(segments), function(i) {
    path <- ancestry(i)
    hits <- path[reset_shift[path]]
    if (!length(hits)) return(NA_character_)
    paste(passage_ids[[hits[1L]]], collapse = '|')
  }, character(1))

  per_segment[[j]] <- data.frame(
    origin = b$origin, key = b$key, cohort = edges$cohort, segment_id = edges$child,
    parent_segment_id = edges$parent, oxygen_pct = edges$oxygen_pct,
    passage_ids = vapply(passage_ids, paste, character(1), collapse = '|'),
    count_target = edges$target_effective, target_source = edges$target_source,
    reference_score_day = base$selection$day, fitted_score_day = fitted$selection$day,
    reference_parent_day = base$selection$parent_day, fitted_parent_day = fitted$selection$parent_day,
    reference_reseed_mode = base$reseed_mode, fitted_reseed_mode = fitted$reseed_mode,
    score_day_changed = obs_shift, incoming_parent_day_changed = reset_shift,
    incoming_reseed_mode_changed = mode_shift, ancestral_parent_day_changes = path_reset,
    ancestral_reseed_mode_changes = path_mode, first_changed_parent_day_child = first_reset,
    reference_target_choice_margin_cells = base$choice_margin_cells,
    fitted_target_choice_margin_cells = fitted$choice_margin_cells,
    start_composition_tv = vapply(seq_along(segments), function(i) {
      tv(normalized(fitted$snapshots[[i]][1, ]), normalized(base$snapshots[[i]][1, ]))
    }, numeric(1)))

  growth_rebuild_error <- 0
  per_growth[[j]] <- do.call(rbind, lapply(names(pred$growth), function(id) {
    i <- passage_lookup[[id]]
    ref_index <- base$selection$index[i]
    fit_index <- fitted$selection$index[i]
    g <- unname(pred$growth[[id]])
    g0 <- unname(case$data$growth[[id]])
    reconstructed <- growth_at(fitted, i, fit_index)
    reconstructed0 <- growth_at(base, i, ref_index)
    growth_rebuild_error <<- max(growth_rebuild_error, abs(reconstructed - g), abs(reconstructed0 - g0))
    same_day <- growth_at(fitted, i, ref_index)
    error <- (g - g0) / growth_scale
    data.frame(origin = b$origin, key = b$key, passage_id = id,
               reference_score_day = base$selection$day[i], fitted_score_day = fitted$selection$day[i],
               reference_parent_day = base$selection$parent_day[i], fitted_parent_day = fitted$selection$parent_day[i],
               ancestral_parent_day_changes = path_reset[i], first_changed_parent_day_child = first_reset[i],
               ancestral_reseed_mode_changes = path_mode[i], generating_growth = g0, fitted_growth = g,
               fitted_growth_at_reference_day = same_day, standardized_error = error,
               standardized_error_at_reference_day = (same_day - g0) / growth_scale,
               loss_contribution = error^2 / (2 * length(case$data$growth)))
  }))

  kary_rebuild_error <- 0
  decomposition_error <- 0
  identity_error <- 0
  width <- candidate$params$sigma_kary
  per_kary[[j]] <- do.call(rbind, lapply(names(pred$kary), function(id) {
    i <- passage_lookup[[id]]
    ref_index <- base$selection$index[i]
    fit_index <- fitted$selection$index[i]
    p0 <- reference$latent_kary[[id]]
    p <- pred$latent_kary[[id]]
    q0 <- case$data$prob[[id]]
    q <- pred$kary[[id]]
    q00 <- observe(p0, truth$sigma_kary)
    q11 <- observe(p, width)
    q10 <- observe(p, truth$sigma_kary)
    q01 <- observe(p0, width)
    common_day_p <- normalized(fitted$snapshots[[i]][ref_index, ])
    common_day_q <- observe(common_day_p, width)
    latent_rebuild <- max(tv(p, normalized(fitted$snapshots[[i]][fit_index, ])),
                          tv(p0, normalized(base$snapshots[[i]][ref_index, ])))
    kary_rebuild_error <<- max(kary_rebuild_error, tv(q00, q0), tv(q11, q), latent_rebuild)
    identity_error <<- max(identity_error, tv(q00, p0), tv(q10, p))
    latent_delta <- q10 - q00
    kernel_delta <- q11 - q10
    decomposition_error <<- max(decomposition_error, max(abs((q - q0) - (latent_delta + kernel_delta))))
    observed_tv <- tv(q, q0)
    fixed_kernel_tv <- tv(q10, q0)
    replay_h2 <- hellinger2(q / sum(q), q0)
    data.frame(origin = b$origin, key = b$key, passage_id = id, fit_width = width,
               generating_width = truth$sigma_kary, reference_score_day = base$selection$day[i],
               fitted_score_day = fitted$selection$day[i], reference_parent_day = base$selection$parent_day[i],
               fitted_parent_day = fitted$selection$parent_day[i], ancestral_parent_day_changes = path_reset[i],
               first_changed_parent_day_child = first_reset[i], ancestral_reseed_mode_changes = path_mode[i],
               start_composition_tv = per_segment[[j]]$start_composition_tv[i],
               observed_tv = observed_tv, latent_tv = tv(p, p0),
               fitted_latent_generating_kernel_tv = fixed_kernel_tv,
               generating_latent_fitted_kernel_tv = tv(q01, q0),
               fitted_latent_kernel_change_tv = tv(q11, q10),
               fitted_kernel_tv_reduction = fixed_kernel_tv - observed_tv,
               fitted_latent_error_at_reference_day_tv = tv(common_day_p, p0),
               fitted_observed_error_at_reference_day_tv = tv(common_day_q, q0),
               generating_latent_mass_below_N66 = sum(p0[grid < 66]),
               fitted_latent_mass_below_N66 = sum(p[grid < 66]),
               fitted_observed_mass_below_N66 = sum(q[grid < 66]),
               fitted_observer_hellinger2 = replay_h2,
               generating_kernel_hellinger2 = hellinger2(q10, q0),
               loss_contribution = replay_h2 / length(case$data$prob))
  }))

  g <- per_growth[[j]]
  k <- per_kary[[j]]
  dg <- sum(g$loss_contribution)
  dk <- sum(k$loss_contribution)
  score_error <- max(abs(dg - b$D_G), abs(dk - b$D_K), abs(dg + dk - b$D),
                     abs(max(k$observed_tv) - b$max_sample_karyotype_TV),
                     abs(max(abs(g$standardized_error)) - b$max_growth_error_over_reference_SD),
                     abs(sqrt(mean(g$standardized_error^2)) - b$growth_RMS_over_reference_SD))
  stopifnot(growth_rebuild_error <= 1e-10, kary_rebuild_error <= 1e-12,
            decomposition_error <= 1e-12, identity_error <= 1e-12, score_error <= 1e-10)
  summaries[[j]] <- data.frame(
    origin = b$origin, key = b$key, D = dg + dk, D_G = dg, D_K = dk,
    karyotype_share_of_D = dk / (dg + dk), score_day_changes = sum(obs_shift),
    incoming_parent_day_changes = sum(reset_shift), incoming_reseed_mode_changes = sum(mode_shift),
    scored_growth_rows_with_ancestral_parent_day_change = sum(g$ancestral_parent_day_changes > 0L),
    karyotype_rows_with_ancestral_parent_day_change = sum(k$ancestral_parent_day_changes > 0L),
    worst_growth_passage = g$passage_id[which.max(abs(g$standardized_error))],
    worst_karyotype_passage = k$passage_id[which.max(k$observed_tv)],
    karyotype_rows_above_tv_gate = sum(k$observed_tv > .005),
    fitted_observed_mean_tv = mean(k$observed_tv), fitted_observed_max_tv = max(k$observed_tv),
    generating_kernel_mean_tv = mean(k$fitted_latent_generating_kernel_tv),
    generating_kernel_max_tv = max(k$fitted_latent_generating_kernel_tv),
    fitted_kernel_improves_tv_rows = sum(k$fitted_kernel_tv_reduction > 1e-12),
    D_K_using_generating_kernel = mean(k$generating_kernel_hellinger2),
    D_K_at_reference_score_days = mean(vapply(seq_len(nrow(k)), function(n) {
      id <- k$passage_id[n]
      i <- passage_lookup[[id]]
      same_day <- normalized(fitted$snapshots[[i]][base$selection$index[i], ])
      hellinger2(observe(same_day, width), case$data$prob[[id]])
    }, numeric(1))),
    D_G_at_reference_score_days = mean(g$standardized_error_at_reference_day^2) / 2)
  validation[[j]] <- data.frame(origin = b$origin, status = 'PASS',
    growth_state_rebuild_error = growth_rebuild_error, karyotype_state_and_kernel_rebuild_tv = kary_rebuild_error,
    generating_kernel_vs_identity_tv = identity_error, residual_decomposition_error = decomposition_error,
    score_metric_rebuild_error = score_error, fitted_reseed_composition_rebuild_tv = fitted$reseed_state_error,
    generating_reseed_composition_rebuild_tv = base$reseed_state_error,
    score_days_checked = length(segments), parent_days_checked = sum(!is.na(parent_index)),
    model_propagations = 0L, optimizer_calls = 0L)
}

dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(out))
options(digits = 17)
write_tsv <- function(name, rows) {
  write.table(rows, file.path(out, name), sep = '\t', quote = TRUE, row.names = FALSE, na = 'NA')
}
write_tsv('input_manifest.tsv', index[match(used, index$path), ])
write_tsv('origin_summary.tsv', do.call(rbind, summaries))
write_tsv('passage_branch_comparison.tsv', do.call(rbind, per_segment))
write_tsv('growth_residuals.tsv', do.call(rbind, per_growth))
write_tsv('karyotype_kernel_crosscheck.tsv', do.call(rbind, per_kary))
write_tsv('saved_search_stage_values.tsv', do.call(rbind, search_stages))
write_tsv('validation.tsv', do.call(rbind, validation))
writeLines(c(R.version.string, 'Execution scope: saved numeric state, observer and arithmetic only.',
             'New model propagations: 0', 'New optimizer calls: 0'), file.path(out, 'runtime.txt'))
cat('Archived D19 attribution passed for all six best found and replay-verified endpoints.\n')
