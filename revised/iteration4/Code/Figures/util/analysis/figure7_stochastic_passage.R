# Integer-day threshold-triggered stochastic passage, isolated Figure7 v4.
# External model kinetics are unchanged; this is NOT the old fitting selector.
f7s_config <- function() {
  reps <- as.integer(Sys.getenv("FIGURE7_STOCHASTIC_REPLICATES", "20"))
  stopifnot(length(reps) == 1L, !is.na(reps), reps >= 2L, reps <= 10000L)
  list(master_seed = 20260904L, rng_kind = "L'Ecuyer-CMRG", replicates = reps,
    checkpoint_days = 1000L, mcse_target_N = 0.01,
    sampling = "systematic_randomized_integerization_then_multivariate_hypergeometric",
    day_convention = "first_eligible_integer_day_actual_population_post_passage_mean")
}

f7s_rng_key <- function(family, seed, oxygen, p, initial) paste(
  "invitro", family, seed, sprintf("%.6f", oxygen),
  sprintf("%.6f", p), initial, sep = "|")

f7s_prepare <- function(endpoint_manifest, run_paths) {
  config <- f7s_config()
  expanded <- endpoint_manifest$expanded
  expanded <- expanded[expanded$model_context == "in vitro", ]
  stopifnot(!anyDuplicated(expanded[c("pair_label", "seed_number")]))
  catalog_path <- file.path(run_paths$run_root, "stochastic_rng_stream_catalog.rds")
  if (file.exists(catalog_path)) {
    catalog <- readRDS(catalog_path)
    stopifnot(identical(catalog$master_seed, config$master_seed),
      identical(catalog$endpoint_identity, expanded[c("pair_label", "seed_number", "parameter_signature")]))
  } else {
    rows <- expand.grid(endpoint = seq_len(nrow(expanded)), oxygen = f7g_o2("in vitro"),
      p = f7ft_p_values(), initial = f7ft_initial_ploidy(), KEEP.OUT.ATTRS = FALSE)
    keys <- with(rows, f7s_rng_key(expanded$pair_label[endpoint],
      expanded$seed_number[endpoint], oxygen, p, initial))
    keys <- sort(keys, method = "radix")
    stopifnot(!anyDuplicated(keys))
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    set.seed(config$master_seed)
    stream <- .Random.seed
    states <- matrix(NA_integer_, 7L, length(keys), dimnames = list(NULL, keys))
    for (i in seq_along(keys)) {
      states[, i] <- stream
      stream <- parallel::nextRNGStream(stream)
    }
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    catalog <- list(master_seed = config$master_seed, states = states,
      endpoint_identity = expanded[c("pair_label", "seed_number", "parameter_signature")],
      rule = "sorted condition key owns stream; replicate r owns nextRNGSubStream^(r-1); independent of workers/chunks/run_id")
    f7ft_atomic_save_rds(catalog, catalog_path, compress = "gzip")
  }
  config$stream_catalog_md5 <- unname(tools::md5sum(catalog_path))
  config_path <- file.path(run_paths$run_root, "stochastic_config.rds")
  if (file.exists(config_path)) stopifnot(identical(readRDS(config_path), config))
  else f7ft_atomic_save_rds(config, config_path)
  f7ft_atomic_write_tsv(data.frame(key = names(config), value = unlist(config)),
    file.path(run_paths$run_root, "stochastic_config.tsv"))
  list(config = config, catalog_path = catalog_path, config_path = config_path,
    catalog = catalog)
}

f7s_streams <- function(catalog, family, seeds, oxygen, p, initial, replicates) {
  indices <- match(f7s_rng_key(family, seeds, oxygen, p, initial),
    colnames(catalog$states))
  if (anyNA(indices)) stop("Condition has no registered random stream.")
  out <- matrix(NA_integer_, 7L, length(seeds) * replicates)
  for (j in seq_along(seeds)) {
    state <- catalog$states[, indices[j]]
    for (r in seq_len(replicates)) {
      out[, (j-1L) * replicates+r] <- state
      state <- parallel::nextRNGSubStream(state)
    }
  }
  out
}

f7s_moments <- function(trajectories, multiplicity, replicates) {
  means <- vars <- matrix(0, multiplicity, ncol(trajectories))
  for (j in seq_len(multiplicity)) {
    x <- trajectories[(j-1L)*replicates+seq_len(replicates), , drop = FALSE]
    means[j, ] <- colMeans(x)
    vars[j, ] <- colSums((x-rep(means[j, ], each = replicates))^2)/(replicates-1L)
  }
  list(sum = colSums(means), sum_squared_endpoint_mean = colSums(means^2),
    sum_within_variance = colSums(vars), sum_mc_variance = colSums(vars)/replicates)
}

# Checkpoints are per parameter identity and O2. Once aggregated into its task,
# a checkpoint is removed; task checkpoints retain completed oxygen conditions.
f7s_operator <- function(step, initial, ploidy, endpoint, oxygen, p, days,
    bundle, checkpoint_path = NULL, fingerprint = "test", keep_trace = FALSE,
    replicates = bundle$stochastic$config$replicates) {
  seeds <- as.integer(strsplit(endpoint$represented_seed_numbers[[1]], ",", fixed=TRUE)[[1]])
  weight <- length(seeds)
  stopifnot(weight == endpoint$endpoint_multiplicity_q10[[1]])
  config <- bundle$stochastic$config
  seed <- bundle$rule$seed_cells[[1]]
  target <- bundle$rule$target_live_cells[[1]]
  empty <- matrix(0, ncol(initial), days+1L)
  state <- list(fingerprint = fingerprint, sum = empty,
    sum_squared_endpoint_mean = empty, sum_within_variance = empty,
    sum_mc_variance = empty, initial_index = 1L, day = 0L,
    active = NULL, summary = list(), traces = list())
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    state <- readRDS(checkpoint_path)
    stopifnot(identical(state$fingerprint, fingerprint))
  }
  # Before any bottleneck all replicas are identical. A deterministic scan
  # detects conditions that never passage and avoids redundant random replicas.
  scan <- f7g_propagate_operator_cpp(step, initial, ploidy, days,
    log(seed), log(target), TRUE)
  for (i in seq_len(ncol(initial))) {
    if (i < state$initial_index) next
    rng <- f7s_streams(bundle$stochastic$catalog, endpoint$pair_label[[1]], seeds,
      oxygen, p, f7ft_initial_ploidy()[i], replicates)
    if (is.null(state$active)) state$active <- list(
      composition = initial[, rep(i, weight*replicates), drop=FALSE],
      log_population = rep(log(seed), weight*replicates), rng = rng,
      count = integer(weight*replicates), first = rep(NA_integer_, weight*replicates),
      last = rep(NA_integer_, weight*replicates), jump = numeric(weight*replicates),
      events = list(), trajectories = list())
    if (isTRUE(scan$no_crossing[[i]])) {
      # No stochastic event: exactly the deterministic trajectory for every
      # original endpoint and replicate, including the same unconsumed streams.
      state$sum[i, ] <- weight * scan$mean_ploidy[i, ]
      state$sum_squared_endpoint_mean[i, ] <- weight * scan$mean_ploidy[i, ]^2
      if (keep_trace) state$active$trajectories <- list(data.frame(day=0:days,
        mean=scan$mean_ploidy[i, ], sd=0, replicate_1=scan$mean_ploidy[i, ],
        replicate_2=scan$mean_ploidy[i, ], replicate_3=scan$mean_ploidy[i, ]))
      state$day <- days
    }
    while (state$day < days) {
      block <- min(config$checkpoint_days, days-state$day)
      active <- state$active
      response <- f7s_propagate_cpp(step, active$composition, ploidy, block,
        seed, target, active$rng, active$log_population, state$day, TRUE,
        if (keep_trace) min(3L, ncol(active$composition)) else 0L)
      moments <- f7s_moments(response$mean_ploidy, weight, replicates)
      indices <- state$day + seq_len(block+1L)
      for (field in names(moments)) state[[field]][i, indices] <- moments[[field]]
      active$composition <- response$final_state
      active$log_population <- response$final_log_population
      active$rng <- response$rng_final
      active$count <- active$count + response$passage_count
      active$first[is.na(active$first)] <- response$first_passage_day[is.na(active$first)]
      active$last[!is.na(response$last_passage_day)] <- response$last_passage_day[!is.na(response$last_passage_day)]
      active$jump <- pmax(active$jump, response$maximum_boundary_mean_jump)
      if (keep_trace) {
        active$events[[length(active$events)+1L]] <- response$events
        active$trajectories[[length(active$trajectories)+1L]] <- data.frame(
          day = state$day+seq.int(0, block),
          mean = colMeans(response$mean_ploidy),
          sd = sqrt(colMeans((response$mean_ploidy-rep(colMeans(response$mean_ploidy),
            each=nrow(response$mean_ploidy)))^2)),
          replicate_1 = response$mean_ploidy[1, ],
          replicate_2 = response$mean_ploidy[2, ],
          replicate_3 = response$mean_ploidy[3, ])
      }
      state$active <- active
      state$day <- state$day+block
      if (!is.null(checkpoint_path)) f7ft_atomic_save_rds(state, checkpoint_path, compress=FALSE)
    }
    a <- state$active
    rows <- lapply(seq_along(seeds), function(j) {
      z <- (j-1L)*replicates+seq_len(replicates)
      data.frame(pair_label=endpoint$pair_label[[1]], endpoint_group=endpoint$endpoint_group[[1]],
        optimizer_seed=seeds[j], endpoint_multiplicity_q10=1L, n_replicate=replicates,
        p_misseg=p, O2_pct=oxygen, initial_ploidy=f7ft_initial_ploidy()[i],
        passage_count=mean(a$count[z]), no_crossing=mean(a$count[z]==0L),
        first_passage_day=if (all(is.na(a$first[z]))) NA_real_ else min(a$first[z], na.rm=TRUE),
        last_passage_day=if (all(is.na(a$last[z]))) NA_real_ else max(a$last[z], na.rm=TRUE),
        maximum_boundary_mean_jump=max(a$jump[z]), protocol_failure_day=NA_integer_,
        maximum_pre_post_mean_error=0, selected_model_day_sum=NA_real_,
        earlier_than_segment_end_count=0L, selected_relative_target_distance_sum=NA_real_)
    })
    state$summary[[i]] <- do.call(rbind, rows)
    if (keep_trace) state$traces[[i]] <- list(initial_ploidy=f7ft_initial_ploidy()[i],
      events=as.data.frame(data.table::rbindlist(a$events)),
      trajectories=as.data.frame(data.table::rbindlist(a$trajectories)),
      continuous=scan$mean_ploidy[i, ], rng_initial=rng, rng_final=a$rng)
    state$initial_index <- i+1L; state$day <- 0L; state$active <- NULL
    if (!is.null(checkpoint_path)) f7ft_atomic_save_rds(state, checkpoint_path, compress=FALSE)
  }
  state$summary <- do.call(rbind, state$summary)
  state
}

f7s_compute_task <- function(task, endpoint_manifest, objective_bundle, contexts,
    paths, run_paths, passage_bundle, fingerprint, smoke=FALSE) {
  output <- task$cache_path[[1]]
  if (file.exists(output)) {
    old <- readRDS(output)
    stopifnot(identical(old$fingerprint, fingerprint), isTRUE(old$qc$passed))
    return(old$qc)
  }
  f7r_load_response_engine(paths)
  indices <- as.integer(strsplit(task$endpoint_indices[[1]], ",", fixed=TRUE)[[1]])
  endpoints <- endpoint_manifest$endpoints[match(indices, endpoint_manifest$endpoints$endpoint_index), ]
  represented <- sum(endpoints$endpoint_multiplicity_q10)
  oi <- seq.int(task$o2_index_start[[1]], task$o2_index_end[[1]])
  oxygen <- f7g_o2("in vitro", smoke)[oi]; days <- f7g_days(smoke)
  reuse <- Sys.getenv("FIGURE7_REUSE_CONTINUOUS_RUN", "")
  if (!nzchar(reuse)) stop("v4 requires an explicitly validated continuous baseline.")
  source_path <- file.path(reuse, "task_cache", basename(output))
  source <- readRDS(source_path)
  keys <- setdiff(names(task), "cache_path")
  model_hash <- function(x) strsplit(x, "|", fixed=TRUE)[[1]][2]
  stopifnot(isTRUE(all.equal(task[keys], source$task[keys], check.attributes=FALSE)),
    identical(model_hash(fingerprint), model_hash(source$fingerprint)),
    identical(days, source$day_values), isTRUE(all.equal(oxygen, source$o2_values)),
    all(is.finite(source$weighted_mean)), isTRUE(source$qc$passed))
  old_endpoints <- f7r_read_tsv(file.path(reuse, "q10_unique_endpoint_manifest.tsv"))
  key <- function(x) paste(x$model_context, x$pair_label, x$endpoint_group)
  ix <- match(key(endpoints), key(old_endpoints))
  stopifnot(!anyNA(ix), identical(endpoints$parameter_signature, old_endpoints$parameter_signature[ix]),
    all(endpoints$endpoint_multiplicity_q10 == old_endpoints$endpoint_multiplicity_q10[ix]))
  shape <- c(5L, length(days), length(oxygen))
  checkpoint <- paste0(output, ".checkpoint.rds")
  work <- list(fingerprint=fingerprint, next_o2=1L, next_endpoint=1L,
    sums=array(0, shape), within=array(0, shape), mcvar=array(0, shape),
    endpoint_square=array(0, shape), summaries=list(), formula_error=0)
  if (file.exists(checkpoint)) {
    work <- readRDS(checkpoint)
    stopifnot(identical(work$fingerprint, fingerprint))
  }
  operator_checkpoint <- paste0(output, ".active_operator.rds")
  for (o in seq_along(oxygen)) {
    if (o < work$next_o2) next
    for (e in seq_len(nrow(endpoints))) {
      if (e < work$next_endpoint) next
      endpoint <- endpoints[e, , drop=FALSE]
      prepared <- f7ft_prepare_endpoint(endpoint, objective_bundle, contexts)
      p <- task$p_misseg[[1]]
      forced <- figure7_force_p_misseg(prepared$run_params, p)
      formula <- figure7_p_misseg_formula_qc(prepared$run_params, p)
      work$formula_error <- max(work$formula_error, formula$maximum_direct_formula_error)
      fixed <- fixo2_fixed_matrix(globalenv(), prepared$config, forced, O2=oxygen[o])
      unit <- prepared$config$N_UNIT %||% 22L
      initial <- f7ft_initial_matrix(fixed$ngrid, unit, f7ft_initial_ploidy())
      operator_key <- paste(fingerprint, o, e, sep="|")
      if (file.exists(operator_checkpoint)) {
        prior_key <- readRDS(operator_checkpoint)$fingerprint
        if (!identical(prior_key, operator_key)) {
          stopifnot(identical(prior_key, work$last_committed_operator_key))
          unlink(operator_checkpoint)
        }
      }
      trace <- e == 1L && oxygen[o] %in% c(.5, 20) && p %in% c(.005, .3)
      response <- f7s_operator(as.matrix(Matrix::expm(fixed$M)), initial,
        as.numeric(fixed$ngrid)/unit, endpoint, oxygen[o], p, max(days),
        passage_bundle, operator_checkpoint, operator_key, trace)
      work$sums[,,o] <- work$sums[,,o] + response$sum
      work$within[,,o] <- work$within[,,o] + response$sum_within_variance
      work$mcvar[,,o] <- work$mcvar[,,o] + response$sum_mc_variance
      work$endpoint_square[,,o] <- work$endpoint_square[,,o] + response$sum_squared_endpoint_mean
      work$summaries[[length(work$summaries)+1L]] <- response$summary
      if (trace) f7ft_atomic_save_rds(list(endpoint=endpoint, O2_pct=oxygen[o], p_misseg=p,
        fingerprint=fingerprint, traces=response$traces),
        file.path(run_paths$run_root, sprintf("stochastic_trace_%s_o%03d_e%03d.rds",
          task$task_id[[1]], oi[o], e)), compress="gzip")
      work$next_endpoint <- e+1L
      work$last_committed_operator_key <- operator_key
      f7ft_atomic_save_rds(work, checkpoint, compress=FALSE)
      # This exact task-owned active checkpoint has been atomically reduced.
      if (file.exists(operator_checkpoint)) unlink(operator_checkpoint)
    }
    work$next_o2 <- o+1L; work$next_endpoint <- 1L
    f7ft_atomic_save_rds(work, checkpoint, compress=FALSE)
    message(task$task_id[[1]], " oxygen ", o, "/", length(oxygen), " complete")
  }
  mean <- work$sums/represented
  day0 <- max(abs(mean[,1,]-f7ft_initial_ploidy()))
  mcse <- sqrt(pmax(work$mcvar, 0))/represented
  qc <- source$qc
  qc$cache_path <- output; qc$continuous_source <- source_path
  qc$maximum_day0_abs_error <- day0; qc$maximum_direct_formula_error <- work$formula_error
  qc$maximum_mcse_N <- max(mcse); qc$n_stochastic_replicate <- passage_bundle$stochastic$config$replicates
  qc$all_finite <- all(is.finite(mean)); qc$protocol_mask_valid <- TRUE
  qc$passed <- qc$all_finite && day0 < 1e-10 && work$formula_error < 1e-12
  qc$mcse_target_met <- max(mcse) <= passage_bundle$stochastic$config$mcse_target_N
  passage <- do.call(rbind, work$summaries)
  qc$total_passage_events <- sum(passage$passage_count * passage$n_replicate)
  qc$maximum_passage_pre_post_mean_error <- 0
  f7ft_atomic_save_rds(list(profile=f7g_profile(), fingerprint=fingerprint, task=task,
    oxygen_index=oi, o2_values=oxygen, day_values=days, represented_optimizer_endpoint=represented,
    weighted_mean=source$weighted_mean, passage_weighted_mean=mean,
    passage_feasible_weight=array(as.integer(represented), shape),
    stochastic_within_endpoint_variance=work$within/represented,
    stochastic_mcse=mcse,
    between_endpoint_mean_variance=pmax((work$endpoint_square-represented*mean^2)/(represented-1), 0),
    continuous_source_md5=unname(tools::md5sum(source_path)), passage=passage, qc=qc), output, compress=FALSE)
  if (file.exists(checkpoint)) unlink(checkpoint)
  qc
}
