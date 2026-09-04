# Prefix-nested replicas: R=20/50/100 use exactly the same first R streams.
f7s_pilot <- function(paths, run_paths, endpoint_manifest, objective_bundle,
    contexts, passage_bundle, n_core) {
  endpoints <- endpoint_manifest$endpoints
  rows <- list()
  for (family in f7ft_family_levels()) {
    candidates <- endpoints[endpoints$model_context=="in vitro" & endpoints$pair_label==family,]
    chosen <- candidates[unique(c(1L,ceiling(nrow(candidates)/2))),,drop=FALSE]
    for (e in seq_len(nrow(chosen))) for (oxygen in c(.5,20)) for (p in f7ft_p_values()) {
      rows[[length(rows)+1L]] <- list(endpoint=chosen[e,,drop=FALSE], oxygen=oxygen, p=p)
    }
  }
  for (i in seq_along(rows)) rows[[i]]$case_id <- sprintf("P%03d",i)
  worker <- function(case) {
    started <- proc.time()[[3]]
    e <- case$endpoint
    prepared <- f7ft_prepare_endpoint(e,objective_bundle,contexts)
    fixed <- fixo2_fixed_matrix(globalenv(),prepared$config,
      figure7_force_p_misseg(prepared$run_params,case$p), O2=case$oxygen)
    unit <- prepared$config$N_UNIT %||% 22L
    initial <- f7ft_initial_matrix(fixed$ngrid,unit,f7ft_initial_ploidy())
    step <- as.matrix(Matrix::expm(fixed$M)); ploidy <- as.numeric(fixed$ngrid)/unit
    seed <- passage_bundle$rule$seed_cells[[1]]; target <- passage_bundle$rule$target_live_cells[[1]]
    days <- 10000L
    scan <- f7g_propagate_operator_cpp(step,initial,ploidy,days,log(seed),log(target),TRUE)
    measures <- list()
    for (i in c(1L,3L,5L)) {
      if (scan$no_crossing[i]) {
        trajectory <- matrix(rep(scan$mean_ploidy[i,],each=100L),100L)
        event_count <- 0
      } else {
        streams <- f7s_streams(passage_bundle$stochastic$catalog,e$pair_label[[1]],
          e$representative_seed_number[[1]],case$oxygen,case$p,f7ft_initial_ploidy()[i],100L)
        response <- f7s_propagate_cpp(step,initial[,rep(i,100L),drop=FALSE],ploidy,
          days,seed,target,streams,rep(log(seed),100L),0L,TRUE,0L)
        trajectory <- response$mean_ploidy; event_count <- sum(response$passage_count)
      }
      reference_mean <- colMeans(trajectory)
      for (r in c(20L,50L,100L)) {
        x <- trajectory[seq_len(r),,drop=FALSE]
        mu <- colMeans(x); variance <- colSums((x-rep(mu,each=r))^2)/(r-1L)
        measures[[length(measures)+1L]] <- data.frame(pair_label=e$pair_label[[1]],
          endpoint_group=e$endpoint_group[[1]], optimizer_seed=e$representative_seed_number[[1]],
          O2_pct=case$oxygen,p_misseg=case$p,initial_ploidy=f7ft_initial_ploidy()[i],
          replicates=r, maximum_endpoint_mcse_N=max(sqrt(variance/r)),
          maximum_mean_delta_vs_100_N=max(abs(mu-reference_mean)),
          maximum_within_sd_N=max(sqrt(variance)), passage_events_100=event_count,
          no_crossing=scan$no_crossing[i], n_state=nrow(step), n_day=days)
      }
    }
    out <- do.call(rbind,measures)
    out$case_elapsed_seconds <- proc.time()[[3]]-started
    out$case_id <- case$case_id
    f7ft_atomic_write_tsv(out,file.path(run_paths$run_root,paste0("stochastic_pilot_",case$case_id,".tsv")))
    message(case$case_id," pilot complete; seconds=",round(out$case_elapsed_seconds[1],1))
    out
  }
  message("Stochastic pilot: ",length(rows)," model conditions, 20/50/100 nested repeats; full 10000 days")
  results <- f7ft_parallel_lapply(rows,worker,n_core=n_core)
  results <- do.call(rbind,results)
  f7ft_atomic_write_tsv(results,file.path(run_paths$run_root,"stochastic_pilot_results.tsv"))
  # A stricter per-endpoint check protects the final 50-endpoint mean. Record
  # unresolved variability rather than pretending a finite pilot proves all cells.
  summary <- do.call(rbind,lapply(c(20L,50L,100L),function(r) {
    z <- results[results$replicates==r,]
    data.frame(replicates=r, maximum_endpoint_mcse_N=max(z$maximum_endpoint_mcse_N),
      maximum_mean_delta_vs_100_N=max(z$maximum_mean_delta_vs_100_N),
      passed=all(z$maximum_endpoint_mcse_N<=passage_bundle$stochastic$config$mcse_target_N))
  }))
  f7ft_atomic_write_tsv(summary,file.path(run_paths$run_root,"stochastic_pilot_selection.tsv"))
  eligible <- summary$replicates[summary$passed]
  if (!length(eligible)) stop("Pilot needs more than 100 repeats; inspect convergence before a full run.")
  selected <- min(eligible)
  # Preserve representative full-range trajectories/variance in the final run;
  # runtime estimates are empirical, not an assertion of guaranteed wall time.
  unique_time <- unique(results[c("pair_label","endpoint_group","O2_pct","p_misseg","case_elapsed_seconds")])
  budget <- data.frame(selected_replicates=selected, pilot_cases=nrow(unique_time),
    median_case_seconds_100_replicates=median(unique_time$case_elapsed_seconds),
    maximum_case_seconds_100_replicates=max(unique_time$case_elapsed_seconds),
    pilot_total_worker_seconds=sum(unique_time$case_elapsed_seconds),
    mcse_target_N=passage_bundle$stochastic$config$mcse_target_N, passed=TRUE)
  f7ft_atomic_write_tsv(budget,file.path(run_paths$run_root,"stochastic_pilot_budget.tsv"))
  message("Selected fixed stochastic repeat count: ",selected)
  invisible(budget)
}
