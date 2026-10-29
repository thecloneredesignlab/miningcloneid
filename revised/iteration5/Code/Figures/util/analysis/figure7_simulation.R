# Figure7: own total trajectories -> two passive origin labels -> all products.
source(file.path(.fp_module_directory,"figure7_source_accounting.R"),local=globalenv())
f7_tasks <- function(a) {
  cfg <- a$cfg;parts <- list()
  for(ctx in c("in vivo","in vitro")) {
    oo <- sort(unique(c(if(ctx=="in vivo")cfg$matrix_oxygen_invivo else cfg$matrix_oxygen_invitro,cfg$main_oxygen)))
    parts[[length(parts)+1L]] <- data.table::CJ(pair_label=c("C01","C02"),model_context=ctx,
      O2_pct=oo,p_misseg=cfg$matrix_p_misseg,threshold_N=cfg$matrix_thresholds_N,
      initial_ploidy=cfg$initial_ploidy)
  }
  parts[[length(parts)+1L]] <- data.table::CJ(pair_label=c("C01","C02"),model_context=c("in vivo","in vitro"),
    O2_pct=cfg$main_oxygen,p_misseg=cfg$supp5_p_misseg,threshold_N=cfg$threshold_N,initial_ploidy=cfg$initial_ploidy)
  parts[[length(parts)+1L]] <- data.table::CJ(pair_label=c("C01","C02"),model_context=c("in vivo","in vitro"),
    O2_pct=c(0,5),p_misseg=.3,threshold_N=66L,initial_ploidy=4)
  z <- unique(data.table::rbindlist(parts))
  if(a$profile=="representative")z <- z[threshold_N==66L&
    (O2_pct==.5|(p_misseg==.3&O2_pct%in%c(0,5)&initial_ploidy==4)|
      (model_context=="in vitro"&O2_pct==20&initial_ploidy==4))]
  data.table::setorder(z,model_context,pair_label,O2_pct,p_misseg,threshold_N,initial_ploidy)
  z[,`:=`(task_id=seq_len(.N),condition_mode="fixed_p_misseg")];z
}
f7_direct_operator <- function(l,c,culture) {
  n <- length(l$ngrid);P <- Q <- PU <- matrix(0,0,0);survive <- numeric();error <- 0
  if(!culture) {
    P <- as.matrix(Matrix::expm(l$M));PU <- as.matrix(Matrix::expm(l$M-c$KH))
    aug <- matrix(0,2*n,2*n);aug[seq_len(n),seq_len(n)] <- l$M
    aug[n+seq_len(n),seq_len(n)] <- c$KH
    aug[n+seq_len(n),n+seq_len(n)] <- -diag(l$lambda+l$mu)
    E <- as.matrix(Matrix::expm(aug));Q <- E[n+seq_len(n),seq_len(n),drop=FALSE]
    survive <- exp(-l$lambda-l$mu)
    error <- max(abs(E[seq_len(n),seq_len(n)]-P))
    if(error>1e-9)stop("Augmented direct-origin propagation altered X")
  }
  rows <- matrix(0,2,n)
  if(!culture) {
    aug <- matrix(0,n+2L,n+2L);aug[seq_len(n),seq_len(n)] <- l$M
    aug[n+1L,seq_len(n)] <- c$hl;aug[n+2L,seq_len(n)] <- c$ll
    rows <- as.matrix(Matrix::expm(aug))[n+1:2,seq_len(n),drop=FALSE]
    if(min(rows)< -1e-12)stop("Negative integrated daughter counter operator")
    rows[rows<0] <- 0
  }
  list(op=list(G=l$G,KH=c$KH,P=P,PU=PU,Q=Q,lambda=l$lambda,mu=l$mu,survive=survive),
    count_rows=rows,error=error)
}
f7_diagnostics <- function(l,r,t) {
  low <- l$ngrid<66L;cf <- if(t$model_context=="in vitro")r$final_crowding_factor else rep(1,ncol(r$final_state))
  X <- rowMeans(sweep(r$final_state,2,cf,"*"));U <- rowMeans(sweep(r$final_untagged_state,2,cf,"*"));T <- X-U
  ordinary <- sweep(l$B,2,l$lambda*(1-l$p_wgd),"*");wgd <- sweep(l$W,2,l$lambda*l$p_wgd,"*")
  KH <- ordinary;KH[!low,] <- 0;KH[,low] <- 0
  bins <- floor(l$ngrid/l$N_unit/.25)*.25;heat <- history <- budget <- list()
  for(branch in c("Ordinary division","WGD")) {
    J <- if(branch=="Ordinary division")ordinary else wgd
    first <- if(branch=="Ordinary division")sweep(KH,2,U,"*") else matrix(0,nrow(J),ncol(J))
    blocks <- list(All=sweep(J,2,X,"*"),`Historical H-to-L`=sweep(J,2,T,"*")+first,
      `Never H-to-L`=sweep(J,2,U,"*")-first)
    if(max(abs(blocks[[1]]-blocks[[2]]-blocks[[3]]))>1e-12)stop("Historical diagnostic births do not add up")
    for(group in names(blocks)) {
      h <- data.table::data.table(expand.grid(daughter=bins,source=bins),flux=1000*as.numeric(blocks[[group]]))
      h <- h[,.(flux=sum(flux)),by=.(source,daughter)]
      h <- cbind(t,branch=branch,h)
      if(group=="All")heat[[length(heat)+1L]] <- h else history[[length(history)+1L]] <- cbind(h,lineage=group)
    }
  }
  J <- ordinary+wgd;HL <- colSums(J[low,,drop=FALSE])*(!low);LL <- colSums(J[low,,drop=FALSE])*low
  lowfrac <- mean(colSums(r$final_state[low,,drop=FALSE]))
  for(group in c("All","Historical H-to-L","Never H-to-L")) {
    weight <- switch(group,All=X,`Historical H-to-L`=T,`Never H-to-L`=U)
    unweighted <- switch(group,All=r$final_state,`Historical H-to-L`=r$final_state-r$final_untagged_state,
      `Never H-to-L`=r$final_untagged_state)
    H <- if(group=="Never H-to-L")0 else sum(HL*X)
    L <- sum(LL*weight);remove <- sum(l$lambda[low]*weight[low])
    death <- mean(colSums(sweep(unweighted,1,l$mu,"*")[low,,drop=FALSE]))
    budget[[length(budget)+1L]] <- cbind(t,lineage=group,data.table::data.table(external_HL=H,
      internal_LL=L,maternal_removal=remove,direct_death=death,low_net_change=H+L-remove-death,low_fraction=lowfrac))
  }
  budget <- data.table::rbindlist(budget)
  list(heat=data.table::rbindlist(heat),history_heat=data.table::rbindlist(history),
    budget=budget[lineage=="All"][,lineage:=NULL],history_budget=budget[lineage!="All"])
}
f7_simulate_task <- function(a,b,t) {
  e <- b$endpoints[pair_label==t$pair_label&model_context==t$model_context]
  stopifnot(nrow(e)==1L);pr <- b$prepared[[e$prepared_index]]
  base <- fp_ledger(pr,t$O2_pct,t$p_misseg);culture <- t$model_context=="in vitro";outputs <- list()
  for(beta in b$config$intervention_q) {
    l <- fp_redirect(base,t$threshold_N,beta);c <- fp_birth_coefficients(l,t$threshold_N)
    x0 <- fp_pure_initial(l,t$initial_ploidy) # Intentional Figure7 diagnostic design, including culture.
    p <- f7_direct_operator(l,c,culture)
    if(culture)r <- fp_culture(b,pr,l,x0,e$pair_label,e$seed_number,t$O2_pct,t$p_misseg,
      Q=t$initial_ploidy,threshold=t$threshold_N,tagged=TRUE,daily=TRUE,events=TRUE) else {
      streams <- fp_streams(e$pair_label,e$seed_number,t$O2_pct,t$p_misseg,t$initial_ploidy,1L)
      r <- f7_continuous_cpp(p$op$P,p$op$PU,matrix(x0,ncol=1L),l$ngrid/l$N_unit,colSums(l$M),
        c$hl,c$ll,b$config$day,b$rule$seed_cells[1],b$rule$target_live_cells[1],
        streams$state,streams$label,FALSE,TRUE)
      r$executed_repeats <- 1L;r$negative_state_clipping_count <- 0L
    }
    dt <- if(culture)pr$config$DT else 1;mode <- if(identical(pr$config$crowding,"logistic"))0L else 1L
    h <- f7_history_cpp(p$op,r,x0,c$low,l$ngrid/l$N_unit,c$hl,c$ll,culture,dt,
      isTRUE(pr$config$Crowding),mode,pr$config$K,b$rule$seed_cells[1],b$config$day)
    d <- f7_direct_cpp(p$op,r,x0,c$low,l$ngrid/l$N_unit,c$hl,c$ll,culture,dt,
      isTRUE(pr$config$Crowding),mode,pr$config$K,b$rule$seed_cells[1],b$config$day)
    profiles <- t$threshold_N==b$config$threshold_N&&t$O2_pct%in%b$config$main_oxygen
    daily <- f7_daily_cpp(p$op,r,x0,c$low,l$ngrid/l$N_unit,c$hl,c$ll,p$count_rows,
      culture,dt,isTRUE(pr$config$Crowding),mode,pr$config$K,b$rule$seed_cells[1],
      b$config$day,profiles,as.integer(b$config$snapshot_days))
    meta <- cbind(t,seed_number=e$seed_number,q=beta)
    nominal <- if(culture)b$replicates else 1L
    series <- data.table::rbindlist(list(
      cbind(meta,analysis="instantaneous",metric="flux_percent",source_accounting="daily_viable_low_daughters",
        fp_summaries(daily$high_parent_percent,nominal)),
      cbind(meta,analysis="instantaneous",metric="growth_percent",source_accounting="daily_viable_low_daughters",
        fp_summaries(daily$low_parent_percent,nominal)),
      cbind(meta,analysis="history",metric="flux_percent",fp_summaries(h$lineage_cell_share,nominal)),
      cbind(meta,analysis="history",metric="growth_percent",fp_summaries(h$lineage_growth_share,nominal))),fill=TRUE)
    source_counts_summary <- f7_count_summary(daily,meta,nominal)
    population_summary <- if(profiles)f7_population_summary(daily,meta,nominal) else NULL
    distribution_summary <- if(profiles)f7_distribution_summary(daily,meta,l,c) else NULL
    X <- r$final_state;D <- d$final_direct_state
    all_births <- colSums(sweep(X,1,c$ll,"*"));direct_births <- colSums(sweep(D,1,c$ll,"*"))
    mass <- colSums(X[c$low,,drop=FALSE]);cf <- if(culture)r$final_crowding_factor else rep(1,ncol(X))
    valid <- if(culture)cf*all_births>1e-14 else mass>1e-14&all_births/mass>1e-8
    expected <- ifelse(valid,100*direct_births/all_births,NA_real_)
    got <- d$growth_percent[,b$config$day+1L]
    endpoint_error <- if(any(valid))max(abs(expected[valid]-got[valid])) else 0
    qc <- cbind(meta,data.table::data.table(direct_flux_error=d$reference_flux_error,
      direct_ploidy_error=d$reference_ploidy_error,passage_error=max(d$saved_passage_error,h$maximum_pre_passage_error),
      final_state_error=max(d$final_state_error,h$maximum_final_state_error),direct_subset_error=d$subset_error,
      direct_growth_balance_error=d$growth_balance_error,history_metric_error=h$maximum_reference_metric_error,
      history_balance_error=h$maximum_balance_error,endpoint_growth_error=endpoint_error,
      propagator_error=p$error,negative_state_clipping_count=r$negative_state_clipping_count+h$negative_state_clipping_count,
      direct_roundoff_corrections=d$roundoff_corrections,n_passages=nrow(r$passage_events),
      executed_repeats=r$executed_repeats,initial_direct_mass=sum(d$direct_cell_share[,1],na.rm=TRUE)))
    qc[,`:=`(daily_passage_error=daily$passage_error,daily_final_state_error=daily$final_state_error,
      daily_direct_state_error=max(abs(daily$final_direct_state-D)),daily_subset_error=daily$subset_error,
      daily_complement_error=daily$complement_error,
      daily_day_zero_empty=all(is.na(daily$high_parent_percent[,1]))&&all(is.na(daily$low_parent_percent[,1])))]
    qc[,passed:=direct_flux_error<1e-8&direct_ploidy_error<1e-8&passage_error<1e-8&final_state_error<1e-8&
      direct_subset_error<1e-9&direct_growth_balance_error<1e-12&history_metric_error<1e-8&history_balance_error<1e-9&
      endpoint_growth_error<1e-6&propagator_error<1e-9&negative_state_clipping_count==0L&direct_roundoff_corrections==0L&
      daily_passage_error<1e-8&daily_final_state_error<1e-8&daily_direct_state_error<1e-8&
      daily_subset_error<1e-9&daily_complement_error<1e-10&daily_day_zero_empty]
    if(!all(qc$passed)||!identical(is.na(expected),is.na(got)))stop("Figure7 task numerical check failed")
    diag <- if(beta==0&&t$threshold_N==66L&&t$initial_ploidy==4&&t$p_misseg==.3&&t$O2_pct%in%c(0,5))
      f7_diagnostics(l,r,meta) else NULL
    outputs[[length(outputs)+1L]] <- list(metadata=meta,summary=series,qc=qc,diagnostics=diag,
      trajectories=r,direct=d,history=h,daily_counts=daily[c("HL_per_start_live_cell","LL_per_start_live_cell",
        "log_start_live_cells","high_parent_percent","low_parent_percent")],
      source_counts_summary=source_counts_summary,population_summary=population_summary,
      distribution_summary=distribution_summary,intervention_qc=l$intervention_qc)
  }
  list(outputs=outputs)
}
f7_aggregate <- function(a,b,tasks) {
  outputs <- unlist(lapply(tasks,`[[`,"outputs"),recursive=FALSE)
  z <- data.table::rbindlist(lapply(outputs,`[[`,"summary"));qc <- data.table::rbindlist(lapply(outputs,`[[`,"qc"))
  fp_write(z[analysis=="instantaneous"],file.path(a$run,"tables/instant_time_summary.tsv"))
  fp_write(z[analysis=="history"],file.path(a$run,"tables/history_time_summary.tsv"))
  fp_write(qc,file.path(a$run,"provenance/task_numerical_validation.tsv"))
  for(spec in list(c("source_counts_summary","daily_source_counts_summary.tsv"),
    c("population_summary","low_population_time_summary.tsv"),
    c("distribution_summary","low_population_distribution_summary.tsv")))
    fp_write(data.table::rbindlist(lapply(outputs,`[[`,spec[1]),fill=TRUE),file.path(a$run,"tables",spec[2]))
  diagnostics <- Filter(Negate(is.null),lapply(outputs,`[[`,"diagnostics"))
  for(spec in list(c("heat","source_daughter_flux_binned.tsv"),c("budget","winner_low_population_budget.tsv"),
    c("history_heat","history_source_daughter_flux_binned.tsv"),c("history_budget","history_low_population_budget.tsv")))
    fp_write(data.table::rbindlist(lapply(diagnostics,`[[`,spec[1])),file.path(a$run,"tables",spec[2]))
  events <- data.table::rbindlist(lapply(outputs,function(o) {
    e <- o$passage_events;if(is.null(e)||!nrow(e))return(NULL);cbind(o$metadata,e)
  }),fill=TRUE)
  fp_write(events,file.path(a$run,"provenance/passage_events.tsv"))
  invisible(z)
}
f7_validate <- function(a,b,check) {
  for(analysis in c("instant","history")) {
    z <- data.table::fread(file.path(a$run,"tables",paste0(analysis,"_time_summary.tsv")))
    check(paste0(analysis,"_complete_task_beta_metric_day_grid"),nrow(z)==nrow(b$tasks)*length(b$config$intervention_q)*2L*(b$config$day+1L))
    check(paste0(analysis,"_unique_keys"),!anyDuplicated(z[,.(task_id,q,metric,day)]))
    check(paste0(analysis,"_shares_within_0_100"),all(!is.finite(z$median)|(z$median>= -1e-8&z$median<=100+1e-8)))
    check(paste0(analysis,"_quantiles_ordered"),all(is.na(z$median)|(z$q25<=z$median&z$median<=z$q75)))
  }
  qc <- data.table::fread(file.path(a$run,"provenance/task_numerical_validation.tsv"))
  check("total_history_direct_trajectory_and_passage_agree",all(qc$passed))
  check("no_negative_state_clipping",all(qc$negative_state_clipping_count==0L))
  check("direct_label_initially_empty",all(qc$initial_direct_mass==0))
  check("daily_counters_match_original_states_and_passages",all(qc$daily_passage_error<1e-8&
    qc$daily_final_state_error<1e-8&qc$daily_direct_state_error<1e-8))
  check("daily_repeat_source_shares_sum_100",all(qc$daily_complement_error<1e-10))
  check("no_completed_interval_at_day_zero",all(qc$daily_day_zero_empty))
  z <- data.table::fread(file.path(a$run,"tables/instant_time_summary.tsv"))
  check("explicit_daily_production_definition",all(z$source_accounting=="daily_viable_low_daughters"))
  A <- z[metric=="flux_percent"];B <- z[metric=="growth_percent"]
  data.table::setorder(A,task_id,q,day);data.table::setorder(B,task_id,q,day)
  check("paired_source_validity",identical(A$n_valid,B$n_valid)&&identical(is.na(A$median),is.na(B$median)))
  ok <- is.finite(A$median)
  check("daily_medians_complement",all(abs(A$median[ok]+B$median[ok]-100)<1e-9))
  check("daily_iqr_complements_reversed",all(abs(A$q25[ok]+B$q75[ok]-100)<1e-9)&
    all(abs(A$q75[ok]+B$q25[ok]-100)<1e-9))
  check("daily_initial_values_NA",all(is.na(z[day==0]$median)))
  pop <- data.table::fread(file.path(a$run,"tables/low_population_time_summary.tsv"))
  check("five_explicit_origin_groups",setequal(unique(pop$group),c("all_low","direct_first_generation",
    "other_low","historical_lineage","never_HL_lineage")))
  check("empty_group_means_are_NA",all(is.na(pop[n_valid==0]$median)))
  check("profile_means_within_low_compartment",all(pop[metric=="mean_ploidy"&is.finite(median)]$median<3+1e-12))
  cdf <- data.table::fread(file.path(a$run,"tables/low_population_distribution_summary.tsv"))
  check("within_cell_distribution_range",all(!is.finite(cdf$cdf)|(cdf$cdf>=-1e-9&cdf$cdf<=1+1e-9)))
  check("empty_distribution_is_NA",all(is.na(cdf[n_valid_repeats==0]$cdf)))
  check("continuous_tumor_has_no_passage",all(qc[model_context=="in vivo"]$n_passages==0L))
  for(name in c("source_daughter_flux_binned.tsv","winner_low_population_budget.tsv",
    "history_source_daughter_flux_binned.tsv","history_low_population_budget.tsv")) {
    z <- data.table::fread(file.path(a$run,"tables",name))
    check(paste0(name,"_all_families_contexts_oxygen_anchors"),
      nrow(unique(z[,.(pair_label,model_context,O2_pct)]))==8L)
  }
}
