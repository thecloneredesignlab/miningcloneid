main<-function() {
 args<-commandArgs(TRUE);stopifnot(length(args)==3L);s<-normalizePath(args[1]);task<-args[2];mode<-args[3]
 source(file.path(Sys.getenv('D0_RELEASE'),'common.R'),local=.GlobalEnv)
 branch<-substr(mode,nchar(mode),nchar(mode));stopifnot(branch%in%c('A','B'))
 is_check<-startsWith(mode,'check_');out<-file.path(s,if(is_check)'checks'else'fits',task)
 stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);started<-proc.time()[3]
 tryCatch({
 c<-d0_context(s,task,FALSE);controls<-d0_kv(file.path(s,'tables/search_controls.tsv'))
 if(startsWith(mode,'global_')) {
  seed<-as.integer(controls[[paste0('generation_seed_',branch)]]);set.seed(seed)
  d0_rds(.Random.seed,file.path(out,'rng_initial.rds'))
  population<-matrix(runif(256L*14L),nrow=256L,ncol=14L)
  d0_rds(population,file.path(out,'requested_initial_population.rds'))
  d0_write(list(origin=branch,seed=seed,NP=256L,dimensions=14L,initialization='Full-box Uniform; no injected truth or real fit',data_sha256=c$data_sha),file.path(out,'origin.tsv'))
  deadline<-as.numeric(Sys.time())+as.numeric(controls$global_soft_seconds)
  cl<-d0_cluster(s,task,FALSE,deadline,as.integer(controls$global_worker_max_callbacks));on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  native_fn<-function(x).GlobalEnv$d0_ev$fn(x)
  native<-NULL;termination<-'native_return'
  native<-tryCatch(DEoptim::DEoptim(native_fn,rep(0,14),rep(1,14),control=list(NP=256L,
       itermax=as.integer(controls$global_max_generations),trace=50L,cluster=cl,initialpop=population,reltol=1e-10,steptol=1000L)),
       error=function(e){if(grepl('Finite callback/time window ended',conditionMessage(e),fixed=TRUE)){termination<<-'time_or_callback_limit';NULL}else stop(e)})
  counts<-parallel::clusterCall(cl,function()d0_ev$state$n);parallel::stopCluster(cl)
  if(!is.null(native))d0_rds(native,file.path(out,'native_deoptim.rds'))
  d0_rds(.Random.seed,file.path(out,'rng_final.rds'))
  files<-sort(list.files(file.path(out,'initial_inputs'),pattern='[.]rds$',full.names=TRUE));stopifnot(length(files)==256L)
  actual<-do.call(rbind,lapply(files,readRDS));stopifnot(all(dim(actual)==dim(population)))
  keys<-function(x)apply(x,1,function(v)paste(sprintf('%.17g',v),collapse=','))
  stopifnot(identical(sort(keys(actual)),sort(keys(population))))
  d0_rds(actual,file.path(out,'actual_initial_population.rds'));d0_write(list(pass=TRUE,rows=256L,columns=14L,maximum_boundary_projection=0,data_sha256=c$data_sha),file.path(out,'initial_population_identity.tsv'))
  winner<-d0_best(out);d0_rds(winner,file.path(out,'accepted_candidate.rds'))
  d0_write(list(pass=TRUE,origin=branch,method='DEoptim_parallel_global_only',termination=termination,
      evaluations=sum(unlist(counts)),native_evaluations=if(is.null(native))NA_integer_ else native$optim$nfeval,
      objective=winner$loss$objective,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
 } else if(startsWith(mode,'local_')) {
  if(task=='local_A_resume_01') {
   allowance<-d0_kv(file.path(s,'authorizations/local_A_resume_01/approval.tsv'))
   stopifnot(allowance$authorized=='true',allowance$remaining_callbacks=='11181',allowance$remaining_soft_seconds=='6878')
   controls$local_max_callbacks<-as.integer(allowance$remaining_callbacks)
   controls$local_soft_seconds<-as.numeric(allowance$remaining_soft_seconds)
   d0_write(list(previous_task='local_A',previous_job='21902503',previous_callbacks=2L,previous_core_seconds=22L,
      remaining_callbacks=controls$local_max_callbacks,remaining_soft_seconds=controls$local_soft_seconds,
      same_global_origin=TRUE),file.path(out,'resume_binding.tsv'))
  }
  source<-file.path(s,'fits',paste0('global_',branch),'accepted_candidate.rds');initial<-readRDS(source)
  stopifnot(identical(initial$data_sha256,c$data_sha));deadline<-as.numeric(Sys.time())+as.numeric(controls$local_soft_seconds)
  ev<-d0_evaluator(c,out,task,deadline,as.integer(controls$local_max_callbacks),infeasible_infinite=TRUE)
  start<-d0_unit(c,initial$params);q<-ev$fn(start);stopifnot(abs(q-initial$loss$objective)<=1e-8)
  d0_write(list(pass=TRUE,source_task=paste0('global_',branch),first_score=q,saved_source_score=initial$loss$objective,difference=q-initial$loss$objective,projection=FALSE),file.path(out,'start_identity.tsv'))
  sb<-tryCatch(nloptr::nloptr(start,ev$fn,lb=rep(0,14),ub=rep(1,14),opts=d0_local_options(controls,deadline-as.numeric(Sys.time()))),
      d0_budget_stop=function(e)NULL)
  if(!is.null(sb))d0_rds(sb,file.path(out,'native_sbplx.rds'))
  best<-d0_best(out);remaining<-as.integer(controls$local_max_callbacks)-ev$state$n
  pattern<-NULL
  if(remaining>5L&&as.numeric(Sys.time())<deadline)pattern<-tryCatch(o2review_pattern(ev$fn,best$unit,o2review_directions(c$map),remaining,deadline,mesh=c(.02,.005,.001,.0001)),d0_budget_stop=function(e)NULL)
  if(!is.null(pattern))d0_rds(pattern,file.path(out,'native_pattern.rds'))
  winner<-d0_best(out);d0_rds(winner,file.path(out,'accepted_candidate.rds'))
  d0_write(list(pass=TRUE,origin=branch,method='NLOPT_LN_SBPLX_then_bounded_pattern',
       termination=if(ev$state$n>=as.integer(controls$local_max_callbacks))'callback_limit'else if(as.numeric(Sys.time())>=deadline)'time_limit'else'native_return',
       evaluations=ev$state$n,objective=winner$loss$objective,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
 } else if(is_check) {
  candidates<-list();proof<-list()
  for(stage in c('global_','local_')) {
   source_task<-paste0(stage,branch)
   if(stage=='local_'&&file.exists(file.path(s,'tables/local_sources.tsv'))) {
    aliases<-d0_tab(file.path(s,'tables/local_sources.tsv'));a<-aliases$task[aliases$origin==branch]
    if(length(a))source_task<-a
   }
   source<-file.path(s,'fits',source_task,'accepted_candidate.rds');if(!file.exists(source))next
   old<-readRDS(source);stopifnot(identical(old$data_sha256,c$data_sha));runs<-d0_forward(c,old$params,task);pred<-d0_predict(c,runs,old$params);fresh<-d0_loss(c,pred)
   stopifnot(abs(fresh$objective-old$loss$objective)<=1e-8,isTRUE(all.equal(pred,old$predictions,tolerance=1e-12,check.attributes=TRUE)))
   proof[[length(proof)+1L]]<-data.frame(source_task=source_task,saved=old$loss$objective,fresh=fresh$objective,difference=fresh$objective-old$loss$objective,independent=fresh$independent)
   candidates[[length(candidates)+1L]]<-list(params=old$params,predictions=pred,runs=runs,loss=fresh,source_task=source_task)
  }
  stopifnot(length(candidates)>0L);best<-candidates[[which.min(vapply(candidates,function(z)z$loss$objective,numeric(1)))]]
  d0_rds(best,file.path(out,'verified_candidate.rds'));d0_emit(do.call(rbind,proof),file.path(out,'candidate_checks.tsv'))
  # Truth is opened only after score-based estimate selection has finished.
  truth<-readRDS(file.path(s,'truth/truth_parameters.rds'));p<-c$map$table
  p$truth<-vapply(p$parameter,function(n)truth[[n]],numeric(1));p$estimate<-vapply(p$parameter,function(n)best$params[[n]],numeric(1))
  p$natural_error<-p$estimate-p$truth;p$log_ratio<-ifelse(p$truth>0&p$estimate>0,log(p$estimate/p$truth),NA_real_)
  p$scaled_error<-p$natural_error/(p$upper-p$lower);p$estimated<-p$status=='active'
  d0_emit(p,file.path(out,'parameter_recovery.tsv'))
  source(file.path(c$paths$archive,'implementation/diagnostics.R'),local=.GlobalEnv)
  oxygen<-sort(unique(unlist(lapply(c$bundle$adapters,function(a)vapply(a$segments,function(z)z$oxygen_pct,numeric(1))))))
  response<-list()
  for(o in oxygen) {
   tr<-pilotdiag_rates(c,truth,c$bundle$cfg,o)$rates;er<-pilotdiag_rates(c,best$params,c$bundle$cfg,o)$rates
   for(N in c(44,88)) {
    ti<-tr[tr$N==N,];ei<-er[er$N==N,]
    sv<-function(pp)pp$buffer_smax*exp(-pp$buffer_beta*(44/N)^pp$buffer_n_exp)
    tt<-pilotdiag_offspring_identity(N,ti$p_per_chromosome,sv(truth));ee<-pilotdiag_offspring_identity(N,ei$p_per_chromosome,sv(best$params))
    tv<-c(lambda=ti$lambda,mu=ti$mu,p_per_chromosome=ti$p_per_chromosome,survival=sv(truth),effective_error=tt$effective_error_probability_conditional_on_viability,biological_daughter_yield=tt$biological_viable_daughters_per_ordinary_division,retained_loss_rate=ti$ordinary_retained_loss,retained_gain_rate=ti$ordinary_retained_gain)
    ev<-c(lambda=ei$lambda,mu=ei$mu,p_per_chromosome=ei$p_per_chromosome,survival=sv(best$params),effective_error=ee$effective_error_probability_conditional_on_viability,biological_daughter_yield=ee$biological_viable_daughters_per_ordinary_division,retained_loss_rate=ei$ordinary_retained_loss,retained_gain_rate=ei$ordinary_retained_gain)
    response[[length(response)+1L]]<-data.frame(N=N,oxygen_pct=o,function_name=names(tv),truth=unname(tv),estimate=unname(ev),error=unname(ev-tv))
   }
  }
  d0_emit(do.call(rbind,response),file.path(out,'function_recovery.tsv'))
  gate<-best$loss$max_growth_error<=.01&&best$loss$max_kary_tv<=.01&&best$loss$max_flow_tv<=.01
  d0_write(c(best$loss,list(pass=TRUE,origin=branch,prediction_recovery=gate,
       source_task=best$source_task,model_forward_calls=length(candidates),data_sha256=c$data_sha,
       all_parameters_unique=FALSE,wall_seconds=proc.time()[3]-started)),file.path(out,'result.tsv'))
  source(file.path(Sys.getenv('D0_RELEASE'),'plots.R'),local=.GlobalEnv)
  d0_plot_recovery(p,best$loss,out,branch)
 } else stop('Unsupported fit mode')
 d0_guard(s,out);d0_seal(out)
 },error=function(e){d0_write(list(failed=TRUE,message=conditionMessage(e),wall_seconds=proc.time()[3]-started),file.path(out,'FAILED.tsv'));stop(e)})
}
main()
