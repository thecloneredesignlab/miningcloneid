#!/usr/bin/env Rscript
# Independent checks AFTER the guarded representative calculation. Retained
# results may be read here as regression references, never by the cold worker.
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
code <- dirname(script);workspace <- normalizePath(file.path(code,"../.."),mustWork=TRUE)
for(module in c("figure78_pipeline.R","figure78_model.R","figure7_simulation.R","figure8_empirical.R","figure8_simulation.R"))
  source(file.path(code,"util/analysis",module),local=globalenv())
args <- fp_args();args$profile <- "representative";args$workspace <- workspace
if(is.null(args$data_root))stop("Specify the completed guarded --data-root")
audit <- dirname(args$data_root);checks <- list()
check <- function(name,pass,value=NA_real_,tolerance=NA_real_) {
  checks[[length(checks)+1L]] <<- data.table::data.table(test=name,value=value,tolerance=tolerance,passed=isTRUE(pass))
  message(if(isTRUE(pass))"PASS " else "FAIL ",name)
}
delta <- function(x,y) {
  if(length(x)!=length(y)||!identical(is.na(x),is.na(y)))return(Inf)
  if(!any(is.finite(x)))return(0)
  max(abs(x[is.finite(x)]-y[is.finite(x)]))
}
a7 <- fp_init(7L,file.path(code,"run_Figure7.R"),args);b7 <- fp_bundle(a7)
a8 <- fp_init(8L,file.path(code,"run_Figure8.R"),args);b8 <- fp_bundle(a8)
for(a in list(a7,a8)) {
  v <- data.table::fread(file.path(a$run,"provenance/numerical_validation.tsv"))
  check(paste0("Figure",a$figure,"_all_current_numerical_gates"),all(v$passed),nrow(v))
}
access <- data.table::fread(file.path(audit,"cold_start_access_log.tsv"))
check("guarded_cold_start_read_log_present",nrow(access)>0,nrow(access))
check("no_retained_result_or_archive_in_guarded_calculation",
  !any(grepl("/(arxives|archives|frozen_sources)/",access$path))&&
  !any(startsWith(access$path,file.path(workspace,"data/Figures/Figure7/")))&&
  !any(startsWith(access$path,file.path(workspace,"data/Figures/Figure8/"))))
fp_engine(a7)

# A direct-origin cell exits on its own first division, without a daughter
# replenishment term. Death belongs in cohort survival, not the plotted births.
d <- f7_direct_cohort_step_cpp(c(.2,.3),c(0,0),c(.1,.5),c(0,0),1,.05)
check("direct_label_not_inherited_at_first_division",delta(d,c(.199,.2925))<1e-14,delta(d,c(.199,.2925)),1e-14)
d <- f7_direct_cohort_step_cpp(c(.2),c(.4),c(.5),c(.1),.8,.05)
check("direct_cohort_equation_inflow_and_own_loss",abs(d-(.2+.05*(.8*.4-(.8*.5+.1)*.2)))<1e-14)

q7 <- data.table::fread(file.path(a7$run,"provenance/task_numerical_validation.tsv"))
check("actual_stochastic_culture_passages_exercised",sum(q7[model_context=="in vitro"]$n_passages)>0,
  sum(q7[model_context=="in vitro"]$n_passages))
checked_births <- 0L;direct_vs_flux <- 0;repeat_candidate <- NULL
for(id in b7$tasks$task_id) {
  o <- fp_read(fp_task_path(a7,id));t <- o$task
  e <- b7$endpoints[pair_label==t$pair_label&model_context==t$model_context];pr <- b7$prepared[[e$prepared_index]]
  for(v in o$outputs) {
    l <- fp_redirect(fp_ledger(pr,t$O2_pct,t$p_misseg),t$threshold_N,v$metadata$q)
    c <- fp_birth_coefficients(l,t$threshold_N);r <- v$trajectories
    X <- r$final_state;D <- v$direct$final_direct_state;U <- r$final_untagged_state
    all <- colSums(sweep(X,1,c$ll,"*"));direct <- colSums(sweep(D,1,c$ll,"*"))
    cf <- if(t$model_context=="in vitro")r$final_crowding_factor else rep(1,ncol(X))
    err <- max(abs(cf*all-v$direct$all_low_births[,b7$config$day+1L]),
      abs(cf*direct-v$direct$direct_low_births[,b7$config$day+1L]))
    if(err>1e-12)stop("Independent current-low birth numerator/denominator mismatch")
    if(max(abs(D[!c$low,,drop=FALSE]))>1e-12)stop("Direct-origin label moved into H")
    if(max(D-X)>1e-9||min(D)< -1e-12)stop("Direct-origin cells are not a population subset")
    history <- colSums(sweep(X-U,1,c$ll,"*"))
    if(any(direct>all+1e-12)||any(history>all+1e-12))stop("Gross low birth subset exceeds total")
    if(delta(v$direct$growth_percent,r$flux_percent)>1e-3)direct_vs_flux <- direct_vs_flux+1L
    if(is.null(repeat_candidate)&&t$model_context=="in vitro"&&v$metadata$q==0&&nrow(r$passage_events)>0)
      repeat_candidate <- list(t=t,e=e,pr=pr,l=l,r=r)
    checked_births <- checked_births+1L
  }
}
check("independent_direct_and_historical_low_birth_subsets",checked_births==nrow(b7$tasks)*3L,checked_births)
check("supp8_direct_proliferation_is_not_daily_or_endpoint_flux",direct_vs_flux>0,direct_vs_flux)
if(!is.null(repeat_candidate)) {
  v <- repeat_candidate;t <- v$t;e <- v$e
  plain <- fp_culture(b7,v$pr,v$l,fp_pure_initial(v$l,t$initial_ploidy),e$pair_label,e$seed_number,
    t$O2_pct,t$p_misseg,t$initial_ploidy,threshold=t$threshold_N,tagged=FALSE,daily=TRUE,events=TRUE)
  check("adding_history_label_does_not_change_total_state",identical(plain$final_state,v$r$final_state))
  check("adding_history_label_does_not_change_state_RNG",identical(plain$rng_final,v$r$rng_final))
  check("adding_history_label_does_not_change_passage_dates",identical(plain$passage_events$day,v$r$passage_events$day))
}

# Continuous cases have no repeat-size difference and must reproduce the
# preserved accepted definitions. No old result entered the cold calculation.
extra <- list()
for(spec in list(list(context="in vivo",family="C01",threshold=44L),
  list(context="in vitro",family="C02",threshold=88L))) {
  t <- data.table::copy(b7$tasks[model_context==spec$context&pair_label==spec$family&
    O2_pct==.5&p_misseg==.1&initial_ploidy==4][1])
  t[,`:=`(threshold_N=spec$threshold,task_id=max(b7$tasks$task_id)+length(extra)+1L)]
  o <- f7_simulate_task(a7,b7,t)
  check(paste0("full_time_threshold_",spec$threshold,"_",gsub(" ","",spec$context)),
    all(vapply(o$outputs,function(v)all(v$qc$passed),TRUE)))
  extra[[length(extra)+1L]] <- list(task=t,outputs=o$outputs)
}
fp_save(list(purpose="independent full-time threshold spotchecks",parent_fingerprint=b7$fingerprint,
  outputs=extra),file.path(audit,"independent_threshold_spotchecks.rds"))
# Daily source integration intentionally changes instantaneous A/B. Only the
# explicitly unchanged historical definition has a retained-value parity gate.
for(kind in "history") {
  columns <- c("pair_label","model_context","O2_pct","p_misseg","threshold_N","initial_ploidy","q","metric","day","median","q25","q75")
  ref <- data.table::fread(file.path(workspace,"data/Figures/Figure7/tables",paste0(kind,"_time_summary.tsv")),select=columns)
  now <- data.table::fread(file.path(a7$run,"tables",paste0(kind,"_time_summary.tsv")),select=columns)
  joined <- merge(now[model_context=="in vivo"],ref,by=setdiff(columns,c("median","q25","q75")))
  err <- max(delta(joined$median.x,joined$median.y),delta(joined$q25.x,joined$q25.y),delta(joined$q75.x,joined$q75.y))
  check(paste0("Figure7_",kind,"_continuous_reference_parity"),nrow(joined)==nrow(now[model_context=="in vivo"])&&err<1e-6,err,1e-6)
}

# Recompute every saved culture seed mean from the actual repeat values.
n_repeat_checks <- 0L
for(id in b8$tasks$task_id) {
  o <- fp_read(fp_task_path(a8,id));if(o$task$model_context!="in vitro")next
  for(r in o$repeats) {
    row <- o$data[seed_number==r$seed_number&Q==r$Q]
    expected <- colMeans(r$ploidy);ix <- match(row$day,r$days)
    if(delta(row$mean_ploidy,expected[ix])>1e-14)stop("Repeat averaging hierarchy is incorrect")
    n_repeat_checks <- n_repeat_checks+nrow(row)
  }
}
check("all_culture_seed_means_recomputed_from_repeats",n_repeat_checks>0,n_repeat_checks)
endpoints <- data.table::fread(file.path(a8$run,"seed_data/seed_endpoint_predictions.tsv.gz"))
mixes <- endpoints[model_context=="in vivo"&Q%in%1:5]
pure <- endpoints[model_context=="in vivo"&Q%in%c(0,Inf)]
mix_error <- 0
for(i in seq_len(nrow(mixes))) {
  m <- mixes[i];p <- pure[pair_label==m$pair_label&seed_number==m$seed_number&O2_pct==m$O2_pct&p_setting==m$p_setting]
  two <- p[Q==0];four <- p[is.infinite(Q)];stopifnot(nrow(two)==1L,nrow(four)==1L)
  weight <- stats::plogis(two$log_population_factor-four$log_population_factor-log(m$Q))
  mix_error <- max(mix_error,abs(m$mean_ploidy-(weight*two$mean_ploidy+(1-weight)*four$mean_ploidy)))
}
check("in_vivo_mixtures_use_finite_population_mass",mix_error<1e-8,mix_error,1e-8)

fp_engine(a8)
for(family in c("C01","C02")) {
  e <- b8$endpoints[pair_label==family&model_context=="in vitro"&seed_number==1L]
  pr <- b8$prepared[[e$prepared_index]];p <- pr$native_parameters[["p_misseg"]];l <- fp_ledger(pr,20,p)
  x <- fp_initial(pr,l,1,TRUE);two <- fp_initial(pr,l,0,TRUE);four <- fp_initial(pr,l,Inf,TRUE)
  check(paste0(family,"_culture_mixed_fitted_initial_state"),delta(x,(two+four)/2)<1e-14)
  daily <- fp_culture(b8,pr,l,x,family,1L,20,p,1,daily=TRUE,events=TRUE)
  final <- fp_culture(b8,pr,l,x,family,1L,20,p,1,daily=FALSE,events=TRUE)
  short <- fp_culture(b8,pr,l,x,family,1L,20,p,1,daily=TRUE,events=TRUE,days=100L)
  check(paste0(family,"_daily_storage_preserves_terminal_state"),identical(daily$final_state,final$final_state))
  check(paste0(family,"_daily_storage_preserves_passage_dates"),identical(daily$passage_events$day,final$passage_events$day))
  err <- delta(daily$mean_ploidy[seq_len(short$executed_repeats),101L],short$mean_ploidy[,101L])
  check(paste0(family,"_day100_is_same_trajectory_prefix"),err==0,err,0)
  check(paste0(family,"_threshold_trigger_and_actual_inoculum"),nrow(daily$passage_events)>0&&
    all(daily$passage_events$cells_before>=b8$rule$target_live_cells[1])&&
    all(daily$passage_events$cells_after==b8$rule$seed_cells[1]))
  values <- b8$prepared[[e$prepared_index]]$native_parameters
  expected <- pmin(1,pmax(0,values[["p_mis_base"]]+p*.mu_eff_of_O2(20,pr$run_params,N=l$ngrid)/
    (.mu_eff_of_O2(20,pr$run_params,N=l$ngrid)+values[["k_o_mis"]])))
  check(paste0(family,"_native_vs_effective_missegregation_distinct"),delta(expected,l$p_mis_effective)<1e-12)
}
toy <- data.table::data.table(group="pooled",mean_ploidy=c(0,1,2,100),pair_label=c("C01","C01","C02","C02"),seed_number=c(1,2,1,2))
s <- fp_quantiles(toy,"group")
check("pool_seed_values_not_family_quantiles",s$median==1.5&&s$q25==.75&&s$q75==26.5)

for(ctx in c("invivo","invitro")) {
  ref <- data.table::fread(file.path(workspace,"data/Figures/Figure8/seed_data",paste0(ctx,"_seed_snapshots.tsv.gz")))
  now <- data.table::fread(file.path(a8$run,"seed_data",paste0(ctx,"_seed_snapshots.tsv.gz")))
  joined <- merge(now,ref,by=c("pair_label","seed_number","O2_pct","Q","day"))
  if(ctx=="invivo") {
    err <- delta(joined$mean_ploidy.x,joined$mean_ploidy.y)
    check("Figure8_invivo_deterministic_reference_parity",nrow(joined)>0&&err<1e-8,err,1e-8)
  } else {
    # The no-crossing bound is conservative; these representative culture
    # cases actually passage. Compare matched 200-repeat means, not a 2-repeat
    # mean to a 200-repeat mean, and never pass an empty comparison subset.
    spec <- data.table::data.table(pair_label=c("C01","C02"),seed_number=1L,O2_pct=c(.5,20),Q=c(0,1))
    key <- digest::digest(list(bundle=b8$fingerprint,spec=spec,repeats=200L),algo="sha256")
    cache_path <- file.path(audit,"independent_stochastic_reference.rds")
    cache <- if(file.exists(cache_path))fp_read(cache_path) else NULL
    if(is.null(cache)||!identical(cache$key,key)) {
      results <- list()
      for(i in seq_len(nrow(spec))) {
        s <- spec[i];e <- b8$endpoints[pair_label==s$pair_label&model_context=="in vitro"&seed_number==s$seed_number]
        pr <- b8$prepared[[e$prepared_index]];p <- pr$native_parameters[["p_misseg"]];l <- fp_ledger(pr,s$O2_pct,p)
        started <- proc.time()[3]
        o <- fp_culture(b8,pr,l,fp_initial(pr,l,s$Q,TRUE),s$pair_label,s$seed_number,s$O2_pct,p,s$Q,
          repeats=200L,daily=FALSE,events=FALSE)
        results[[i]] <- list(spec=s,complete_input_signature=e$complete_input_signature,
          repeat_ploidy=as.numeric(o$mean_ploidy[,1]),passage_count=o$passage_count,
          final_state=o$final_state,rng_final=o$rng_final,seconds=unname(proc.time()[3]-started))
      }
      cache <- list(purpose="independent matched-200-repeat regression only",key=key,
        parent_fingerprint=b8$fingerprint,spec=spec,repeats=200L,results=results)
      fp_save(cache,cache_path)
    }
    comparison <- list()
    for(r in cache$results) {
      s <- r$spec;e <- b8$endpoints[pair_label==s$pair_label&model_context=="in vitro"&seed_number==s$seed_number]
      old <- ref[pair_label==s$pair_label&seed_number==s$seed_number&O2_pct==s$O2_pct&Q==s$Q&day==1000]
      err <- delta(mean(r$repeat_ploidy),old$mean_ploidy)
      check(paste0("Figure8_culture_matched_200_repeat_reference_",s$pair_label),
        length(r$repeat_ploidy)==200L&&identical(r$complete_input_signature,e$complete_input_signature)&&
        sum(r$passage_count)>0&&nrow(old)==1L&&err<1e-10,err,1e-10)
      comparison[[length(comparison)+1L]] <- cbind(s,repeats=length(r$repeat_ploidy),
        current_mean=mean(r$repeat_ploidy),reference_mean=old$mean_ploidy,absolute_error=err,
        total_passages=sum(r$passage_count))
    }
    fp_write(data.table::rbindlist(comparison),file.path(audit,"matched_200_repeat_regression.tsv"))
  }
}
for(name in c("empirical_tissue_ploidy_summary.tsv","author_figure5_points.tsv")) {
  ref <- data.table::fread(file.path(workspace,"data/Figures/Figure8/tables",name))
  now <- data.table::fread(file.path(a8$run,"tables",name))
  keys <- if(grepl("empirical",name))c("site","cohort") else "row_id"
  joined <- merge(now,ref,by=keys)
  cols <- if(grepl("empirical",name))c("n","mean","median","q25","q75") else c("wgd_fraction","oxygen","oxygen_sd","oxygen_low","oxygen_high")
  err <- max(vapply(cols,function(n)delta(joined[[paste0(n,".x")]],joined[[paste0(n,".y")]]),0.0))
  check(paste0("raw_author_regeneration_",name),nrow(joined)==nrow(ref)&&err<1e-12,err,1e-12)
}
# Exercise every plotting slot using the newly generated layout. This tests the
# current resolver, not a full-size publication claim for representative data.
source(file.path(code,"util/analysis/figure78_render_only.R"),local=globalenv())
ctx <- new.env(parent=emptyenv());ctx$workspace <- workspace;ctx$data_root <- normalizePath(args$data_root)
ctx$registry <- f78ro_registry();ctx$products <- ctx$registry$stem
ctx$current <- list(Figure7=list(),Figure8=list())
ctx$checked <- data.frame(path=character(),relative_path=character(),sha256=character())
ctx$inputs <- data.table::rbindlist(lapply(c(7L,8L),function(figure) {
  m <- data.table::fread(file.path(ctx$data_root,paste0("Figure",figure),"provenance/plot_inputs.tsv"))
  m[,`:=`(figure=paste0("Figure",figure),group="current")];m
}))
ctx$contexts <- f78ro_contexts(ctx)
for(i in seq_len(nrow(ctx$registry))) {
  row <- ctx$registry[i,]
  paths <- f78ro_data(ctx,row$key,f78ro_required_tables(row))
  check(paste0("newly_generated_plot_slots_",row$stem),all(file.exists(paths)),length(paths))
}
rejected <- tryCatch({f78ro_current_inputs(ctx);FALSE},error=function(e)
  grepl("Only validated full-ensemble",conditionMessage(e),fixed=TRUE))
check("representative_data_cannot_masquerade_as_production",rejected)
for(figure in c(7L,8L)) {
  full <- if(figure==7L)a7 else a8;full$profile <- "full"
  tasks <- if(figure==7L)f7_tasks(full) else f8_tasks(full)
  check(paste0("Figure",figure,"_full_manifest_available_without_simulation"),
    nrow(tasks)==if(figure==7L)320L else 9160L,nrow(tasks))
}
protection_path <- file.path(audit,"protected_before.tsv")
if(!file.exists(protection_path))protection_path <- file.path(dirname(audit),"protected_before.tsv")
protected <- data.table::fread(protection_path)
observed <- unname(vapply(protected$path,digest::digest,"",file=TRUE,algo="sha256"))
check("all_formal_scientific_inputs_and_figures_unchanged",identical(observed,protected$sha256_before),nrow(protected))
out <- data.table::rbindlist(checks)
fp_write(out,file.path(audit,"independent_acceptance.tsv"))
fp_write(data.frame(path=script,sha256=fp_sha(script)),file.path(audit,"independent_acceptance_code.tsv"))
if(!all(out$passed))stop("Independent acceptance failed: ",paste(out[passed==FALSE]$test,collapse=", "))
message("PASS: ",nrow(out)," independent full-pipeline acceptance checks")
