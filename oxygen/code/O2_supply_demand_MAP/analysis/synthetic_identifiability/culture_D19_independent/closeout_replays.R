# Offline audit of sealed records. This script never loads the simulator or optimizer.
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
emit <- function(data, path) {
  old <- options(digits = 17)
  on.exit(options(old))
  write.table(data, path, sep = '\t', quote = TRUE, row.names = FALSE)
}

d19_metrics <- function(case, pred) {
  stopifnot(identical(names(pred$growth), names(case$data$growth)),
    identical(names(pred$kary), names(case$data$prob)), length(pred$growth) == 114L,
    length(pred$kary) == 12L, all(lengths(pred$kary) == 133L),
    all(is.finite(pred$growth)), is.finite(case$known_sigma_growth), case$known_sigma_growth > 0)
  residual <- (pred$growth - case$data$growth) / case$known_sigma_growth
  q <- lapply(pred$kary, function(p) {
    stopifnot(all(is.finite(p)), all(p >= 0), abs(sum(p) - 1) <= 1e-12)
    p / sum(p)
  })
  target <- case$data$prob
  tv <- vapply(names(target), function(id) sum(abs(q[[id]] - target[[id]])) / 2, numeric(1))
  hellinger <- vapply(names(target), function(id) {
    stopifnot(length(target[[id]]) == 133L, abs(sum(target[[id]]) - 1) <= 1e-12)
    sum((sqrt(q[[id]]) - sqrt(target[[id]]))^2) / 2
  }, numeric(1))
  growth <- sum(residual^2) / (2 * length(residual))
  kary <- sum(hellinger) / length(hellinger)
  bounds <- c(D = 1e-4, rms_growth = 0.010, max_growth = 0.030, max_kary_tv = 0.005)
  metrics <- c(D = growth + kary, rms_growth = sqrt(mean(residual^2)),
               max_growth = max(abs(residual)), max_kary_tv = max(tv))
  list(metrics = metrics, components = c(D_G = growth, D_K = kary),
       checks = metrics <= bounds, tv = tv, worst_sample = names(tv)[which.max(tv)])
}

audit_replays <- function(root, out) {
  case_file <- file.path(root, 'inputs/public/D19.rds')
  stopifnot(digest::digest(file = case_file, algo = 'sha256') ==
    '0b1787b85af14c024c420c281c59c50ef92ac5a4f10ada926f1cc614ced1cbfe')
  case <- readRDS(case_file)
  rows <- read_tsv(file.path(root, 'reports/verified_candidates.tsv'))
  jobs <- read_tsv(file.path(root, 'accounting/jobs.tsv'))
  repair_path <- file.path(root, 'accounting/technical_retry.tsv')
  repair <- if (file.exists(repair_path)) read_tsv(repair_path) else NULL
  stopifnot(case$kind == 'D', !anyDuplicated(rows[c('origin','key')]), all(rows$verified))
  summaries <- samples <- prediction <- list()
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, ]
    key <- as.character(row$key)
    stopifnot(grepl('^[0-9a-f]{64}$', key), grepl('^D_[ABDF]$', row$origin),
      startsWith(row$first_process_task, row$origin), startsWith(row$second_process_task, row$origin))
    part <- function(task) file.path('verification', task, paste0('point_', key, '.rds'))
    chosen <- c(part(row$first_process_task), part(row$second_process_task))
    stopifnot(identical(chosen[2L], row$state_file))
    stored <- lapply(chosen, function(path) readRDS(file.path(root, path)))
    pass <- logical(2L)
    job_id <- character(2L)
    for (j in seq_along(stored)) {
      ep <- stored[[j]]
      candidate <- ep$candidate
      sched <- jobs[jobs$task == ep$task & jobs$state == 'COMPLETED' & jobs$exit_code == '0:0', , drop = FALSE]
      stopifnot(nrow(sched) == 1L, ep$process > 0,
        identical(ep$task, c(row$first_process_task, row$second_process_task)[j]),
        identical(candidate$key, key), identical(candidate$origin, row$origin),
        length(candidate$unit) == 19L, length(candidate$params) == 20L,
        identical(candidate$params$sigma_growth, case$known_sigma_growth),
        identical(digest::digest(as.numeric(candidate$unit), algo = 'sha256'), key),
        identical(ep$pred$selection, candidate$pred$selection))
      job_id[j] <- as.character(sched$job_id)
      identity_source <- 'CONTAINER_AND_SCHEDULER'
      if (identical(ep$slurm_job, '')) {
        # Recover the single pre-fix receipt identity from its sealed phase and Slurm ledger.
        stopifnot(j == 1L, ep$task == 'D_B_replay1a', nrow(repair) == 1L,
          repair$reason == 'APPTAINER_CLEANENV_REMOVED_SLURM_JOB_ID',
          as.character(repair$prior_first_job) == job_id[j], repair$resume_from == ep$task)
        identity_source <- 'SEALED_PHASE_AND_SCHEDULER'
      } else {
        stopifnot(identical(ep$slurm_job, job_id[j]))
      }
      a <- d19_metrics(case, ep$pred)
      rep_result <- unlist(row[c(names(a$metrics), names(a$components))], use.names = TRUE)
      expected <- c(a$metrics, a$components)
      stopifnot(all(abs(rep_result - expected) <= 1e-12),
        abs(ep$loss$value - a$metrics[['D']]) <= 1e-12)
      pass[j] <- all(a$checks)
      summaries[[length(summaries) + 1L]] <- data.frame(origin = row$origin, key = key,
        replay = j, task = ep$task, slurm_job = job_id[j], R_process = ep$process,
        identity_source = identity_source,
        as.list(a$metrics), as.list(a$components),
        D_gate = a$checks[['D']], growth_RMS_gate = a$checks[['rms_growth']],
        growth_max_gate = a$checks[['max_growth']], karyotype_gate = a$checks[['max_kary_tv']],
        passes_all_four = pass[j], worst_sample = a$worst_sample)
      samples[[length(samples) + 1L]] <- data.frame(origin = row$origin, key = key, replay = j,
        sample = names(a$tv), TV = unname(a$tv), passes = unname(a$tv <= 0.005))
    }
    x <- stored[[1L]]
    y <- stored[[2L]]
    xpred <- x$pred
    ypred <- y$pred
    stopifnot(job_id[1L] != job_id[2L], x$task != y$task,
      identical(x$candidate, y$candidate), identical(xpred$selection, ypred$selection),
      all(pass == row$prediction_pass), identical(names(xpred$states), names(ypred$states)),
      length(xpred$states) > 0L)
    state_error <- max(vapply(names(xpred$states), function(name) {
      left <- xpred$states[[name]]
      right <- ypred$states[[name]]
      stopifnot(identical(dim(left), dim(right)))
      max(abs(left - right) / pmax(1, abs(left), abs(right)))
    }, numeric(1)))
    growth_error <- max(abs(xpred$growth - ypred$growth))
    kary_error <- max(vapply(names(xpred$kary), function(name) {
      sum(abs(xpred$kary[[name]] - ypred$kary[[name]])) / 2
    }, numeric(1)))
    stopifnot(max(state_error, growth_error, kary_error) <= 1e-10)
    prediction[[i]] <- data.frame(origin = row$origin, key = key, maximum_scaled_state_difference = state_error,
      maximum_growth_difference = growth_error, maximum_karyotype_TV_difference = kary_error,
      identical_selection = TRUE, independent_allocations = TRUE, audit = 'PASS')
  }
  emit(do.call(rbind, summaries), file.path(out, 'd19_replay_gate_audit.tsv'))
  emit(do.call(rbind, samples), file.path(out, 'd19_all_sample_tv_audit.tsv'))
  emit(do.call(rbind, prediction), file.path(out, 'd19_replay_pair_audit.tsv'))
  cat('OFFLINE_REPLAY_AUDIT_PASS\n')
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 2L)
  audit_replays(args[1L], args[2L])
}
