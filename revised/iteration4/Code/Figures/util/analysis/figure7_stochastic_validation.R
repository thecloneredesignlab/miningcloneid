# Independent daily R reference and distribution/reproducibility tests.
f7s_test_streams <- function(n, seed=20260904L) {
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); set.seed(seed)
  result <- matrix(NA_integer_, 7L, n); stream <- .Random.seed
  for (i in seq_len(n)) { result[,i] <- stream; stream <- parallel::nextRNGStream(stream) }
  result
}

f7s_reference <- function(step, initial, grid, days, seed, target, stream) {
  assign(".Random.seed", stream, envir=.GlobalEnv)
  p <- as.numeric(initial); logpop <- log(seed)
  means <- numeric(days+1); means[1] <- sum(p*grid); events <- list()
  for (day in seq_len(days)) {
    x <- as.numeric(step %*% p); growth <- sum(x); p <- x/growth
    logpop <- logpop+log(growth)
    if (logpop >= log(target)) {
      before <- sum(p*grid); actual <- exp(logpop); T <- floor(actual+.5)
      u <- runif(1); cumulative <- cumsum(T*p)
      upper <- pmin(T, floor(cumulative+u)); upper[length(upper)] <- T
      counts <- diff(c(0, upper)); draws <- numeric(length(p)); remaining <- T; needed <- seed
      for (i in seq_along(p)) {
        good <- counts[i]; bad <- remaining-good
        draws[i] <- if (needed<=0 || good<=0) 0 else if (bad<=0) needed else
          if (needed==remaining) good else rhyper(1, good, bad, needed)
        needed <- needed-draws[i]; remaining <- remaining-good
      }
      stopifnot(sum(draws)==seed, all(draws<=counts))
      p <- draws/seed; logpop <- log(seed)
      events[[length(events)+1L]] <- data.frame(day=day, cells_before=actual,
        cells_after=seed, ploidy_before=before, ploidy_after=sum(p*grid))
    }
    means[day+1L] <- sum(p*grid)
  }
  list(mean=means, state=p, events=as.data.frame(data.table::rbindlist(events)),
    rng=.Random.seed, logpop=logpop)
}

f7s_validate <- function(paths, run_paths, actual_cases=list()) {
  records <- list()
  check <- function(name, pass, error=0) {
    records[[length(records)+1L]] <<- data.frame(test=name, error=error, passed=isTRUE(pass))
    if (!isTRUE(pass)) stop("Stochastic validation failed: ", name, " error=", error)
  }
  step <- matrix(c(1.3,.1,.2,1.5), 2, 2)
  initial <- matrix(c(.6,.4), 2, 1); grid <- c(2,4); streams <- f7s_test_streams(6)
  calc <- function(step=step, initial=initial, days=31L, seeds=streams[,1,drop=FALSE],
    size=100, target=251, logpop=rep(log(size), ncol(initial)), start=0L, stochastic=TRUE) {
    f7s_propagate_cpp(step, initial, grid, days, size, target, seeds, logpop,
      start, stochastic, ncol(initial))
  }
  # Explicit arguments avoid self-referential default promises in R.
  run <- calc(step, initial)
  ref <- f7s_reference(step, initial, grid, 31L, 100, 251, streams[,1])
  err <- max(abs(run$mean_ploidy[1,]-ref$mean))
  check("independent_R_daily_reference", err<1e-10, err)
  check("independent_R_event_days", identical(as.integer(run$events$day), as.integer(ref$events$day)))
  check("independent_R_final_rng", identical(as.integer(run$rng_final), ref$rng))
  repeat_run <- calc(step, initial)
  check("exact_repeat", identical(run, repeat_run))
  first <- calc(step, initial, days=13L)
  second <- calc(step, first$final_state, days=18L, seeds=first$rng_final,
    logpop=first$final_log_population, start=13L)
  joined <- cbind(first$mean_ploidy, second$mean_ploidy[,-1,drop=FALSE])
  check("bitwise_checkpoint_continuation", identical(unname(run$mean_ploidy), unname(joined)) &&
    identical(run$rng_final, second$rng_final) && identical(run$final_state, second$final_state))
  batch_initial <- initial[, rep(1,6), drop=FALSE]
  batch <- calc(step, batch_initial, seeds=streams)
  perm <- c(6,1,4,2,5,3)
  reordered <- calc(step, batch_initial[,perm,drop=FALSE], seeds=streams[,perm,drop=FALSE])
  check("worker_order_independent", identical(unname(batch$mean_ploidy),
    unname(reordered$mean_ploidy[order(perm),,drop=FALSE])))
  worker <- function(j) calc(step, initial, seeds=streams[,j,drop=FALSE])$mean_ploidy
  serial <- lapply(1:6, worker)
  parallel <- parallel::mclapply(1:6, worker, mc.cores=2L, mc.set.seed=FALSE)
  check("serial_parallel_equal", identical(serial, parallel))
  doubling <- diag(2,2)
  rounded_day <- calc(doubling, initial, days=5L, size=100, target=501)
  check("ceil_crossing_day_3", rounded_day$first_passage_day==3L &&
    abs(rounded_day$events$cells_before[1]-800)<1e-9)
  check("uses_actual_day_population_not_target", rounded_day$events$cells_before[1]>501)
  slow <- calc(diag(1.02,2), initial, days=100L, size=100, target=251)
  check("waits_beyond_five_days", slow$first_passage_day>5L && all(is.finite(slow$mean_ploidy)))
  nohit <- calc(diag(.99,2), initial, days=100L)
  check("no_crossing_remains_valid", nohit$passage_count==0L && all(is.finite(nohit$mean_ploidy)) &&
    identical(nohit$rng_final, streams[,1,drop=FALSE]))
  scalar <- calc(step, initial, stochastic=FALSE)
  continuous <- f7g_propagate_operator_cpp(step, initial, grid, 31L, log(100), log(251), FALSE)
  err <- max(abs(scalar$mean_ploidy-continuous$mean_ploidy))
  check("scalar_passage_continuous_equivalence", err<1e-12, err)
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260904)
  sample <- replicate(10000, f7s_draw_cpp(c(.25,.75), 1000, 100)$sample[1])
  expected_var <- 100*.25*.75*900/999
  check("hypergeometric_mean", abs(mean(sample)-25)<6*sqrt(expected_var/length(sample)), mean(sample)-25)
  check("hypergeometric_variance", abs(var(sample)/expected_var-1)<.06, var(sample)/expected_var-1)
  rounded <- replicate(10000, f7s_draw_cpp(c(.0004,.9996),1000,100)$population[1])
  check("rare_state_unbiased_integerization", abs(mean(rounded)-.4)<.025, mean(rounded)-.4)
  full <- f7s_draw_cpp(c(.1,.2,.7), 1000, 1000)
  check("full_population_draw_conserves_all_states", identical(full$sample,full$population))
  # The R checkpoint wrapper also restores the original endpoint multiplicity.
  endpoint <- data.frame(pair_label="C01",endpoint_group="synthetic",
    represented_seed_numbers="17,21",endpoint_multiplicity_q10=2L)
  catalog_keys <- as.vector(outer(c(17,21),2:6,function(seed,initial)
    f7s_rng_key("C01",seed,.5,.005,initial)))
  catalog_states <- f7s_test_streams(length(catalog_keys))
  colnames(catalog_states) <- catalog_keys
  config <- f7s_config(); config$replicates <- 3L; config$checkpoint_days <- 3L
  bundle <- list(rule=data.frame(seed_cells=100,target_live_cells=251),
    stochastic=list(config=config,catalog=list(states=catalog_states)))
  step5 <- diag(rep(1.25,5)) + matrix(.03,5,5)
  directory <- file.path(run_paths$run_root,"stochastic_wrapper_tests")
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  checkpoint <- file.path(directory,"operator_checkpoint.rds")
  # Explicit new test-owned filename; no previous scientific result is touched.
  if (file.exists(checkpoint)) unlink(checkpoint)
  wrapped <- f7s_operator(step5,diag(5),2:6,endpoint,.5,.005,17L,bundle,
    checkpoint,"synthetic_v4",TRUE)
  resumed <- f7s_operator(step5,diag(5),2:6,endpoint,.5,.005,17L,bundle,
    checkpoint,"synthetic_v4",TRUE)
  check("completed_operator_resume_equal",identical(wrapped,resumed))
  check("original_duplicate_endpoints_retained", nrow(wrapped$summary)==10L &&
    all(wrapped$summary$n_replicate==3L) && setequal(wrapped$summary$optimizer_seed,c(17,21)))
  check("operator_day_zero", max(abs(wrapped$sum[,1]/2-2:6))<1e-12)
  interrupted_path <- file.path(directory,"interrupted_operator.rds")
  if (file.exists(interrupted_path)) unlink(interrupted_path)
  original_save <- f7ft_atomic_save_rds
  assign("f7ft_atomic_save_rds",function(object,path,...) {
    value <- original_save(object,path,...)
    if (identical(path,interrupted_path) && identical(object$day,3L)) stop("test_checkpoint_interrupt")
    value
  },envir=.GlobalEnv)
  interruption <- try(f7s_operator(step5,diag(5),2:6,endpoint,.5,.005,17L,bundle,
    interrupted_path,"synthetic_v4",TRUE),silent=TRUE)
  assign("f7ft_atomic_save_rds",original_save,envir=.GlobalEnv)
  check("interruption_checkpoint_written",inherits(interruption,"try-error") && file.exists(interrupted_path))
  recovered <- f7s_operator(step5,diag(5),2:6,endpoint,.5,.005,17L,bundle,
    interrupted_path,"synthetic_v4",TRUE)
  check("interrupted_operator_resume_equal",identical(wrapped,recovered))
  for (case in actual_cases) {
    init <- case$initial[,1,drop=FALSE]
    actual <- f7s_propagate_cpp(case$step, init, case$ploidy, case$days,
      round(case$seed), ceiling(case$target), streams[,1,drop=FALSE], log(round(case$seed)), 0L, TRUE, 1L)
    reference <- f7s_reference(case$step, init, case$ploidy, case$days,
      round(case$seed), ceiling(case$target), streams[,1])
    error <- max(abs(actual$mean_ploidy[1,]-reference$mean))
    check(paste0("external_model_daily_reference_",case$name), error<1e-8, error)
  }
  result <- do.call(rbind, records)
  path <- file.path(run_paths$run_root, "stochastic_passage_validation.tsv")
  f7ft_atomic_write_tsv(result,path)
  message("Stochastic passage validation: ",nrow(result)," tests passed")
  invisible(result)
}
