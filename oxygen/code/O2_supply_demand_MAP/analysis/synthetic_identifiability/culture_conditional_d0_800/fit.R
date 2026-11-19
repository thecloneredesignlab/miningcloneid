main<-function() {
 args<-commandArgs(TRUE);stopifnot(length(args)==3L);s<-normalizePath(args[1]);task<-args[2];mode<-args[3]
 source(file.path(Sys.getenv('D0_RELEASE'),'common.R'),local=.GlobalEnv)
 positive<-startsWith(mode,'positive_');branch<-if(positive)'positive'else substr(mode,nchar(mode),nchar(mode));stopifnot(branch%in%c('A','B','C','D','E','F','positive'))
 is_check<-startsWith(mode,'check_')||startsWith(mode,'globalcheck_')||mode=='positive_check';out<-file.path(s,if(is_check)'checks'else'fits',task)
 stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);started<-proc.time()[3]
 tryCatch({
 c<-d0_context(s,task,FALSE);controls<-d0_kv(file.path(s,'tables/search_controls.tsv'))
 if(startsWith(mode,'screen_')) {
  seed<-as.integer(controls[[paste0('generation_seed_',branch)]]);set.seed(seed);d0_rds(.Random.seed,file.path(out,'rng_initial.rds'))
  deadline<-as.numeric(Sys.time())+3300;ev<-d0_evaluator(c,out,task,deadline,2000L,save_native=FALSE)
  draw<-function(i){proposal<-if(runif(1)<.5)'growth_favoring'else'uniform';p<-o2review_draw(c$map,proposal,i);list(unit=d0_unit(c,p),proposal=proposal)}
  draws<-list();proposals<-character();good<-integer();ended<-'eight_feasible_anchors'
  for(i in 1:2000) {
   if(as.numeric(Sys.time())>=deadline){ended<-'screen_time_limit';break}
   z<-draw(i);draws[[i]]<-z$unit;proposals[i]<-z$proposal
   q<-tryCatch(ev$fn(z$unit),d0_budget_stop=function(e)NA_real_)
   if(is.na(q)){ended<-'screen_callback_limit';break}
   if(isTRUE(ev$state$last_feasible))good<-c(good,i)
   if(length(good)>=8L)break
  }
  screened<-length(draws)
  d0_emit(data.frame(draw=seq_len(screened),proposal=proposals,feasible=seq_len(screened)%in%good),file.path(out,'screen_provenance.tsv'))
  if(length(draws)<256L)for(i in (length(draws)+1L):256L){z<-draw(i);draws[[i]]<-z$unit;proposals[i]<-z$proposal}
  bank<-do.call(rbind,draws);d0_rds(bank,file.path(out,'raw_draw_bank.rds'));stopifnot(length(good)>=8L)
  population<-bank[1:256,,drop=FALSE];population[1:8,]<-bank[good[1:8],,drop=FALSE]
  d0_rds(population,file.path(out,'requested_initial_population.rds'));d0_rds(.Random.seed,file.path(out,'rng_final.rds'))
  d0_write(list(pass=TRUE,origin=branch,seed=seed,screened=screened,feasible_anchors=length(good),termination=ended,selection='First feasible in draw order; no objective ranking or truth distances',data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
 } else if(startsWith(mode,'global_')) {
  seed<-as.integer(controls[[paste0('generation_seed_',branch)]]);set.seed(seed)
  d0_rds(.Random.seed,file.path(out,'rng_initial.rds'))
  population<-readRDS(file.path(s,'fits',paste0('screen_',branch),'requested_initial_population.rds'))
  d0_rds(population,file.path(out,'requested_initial_population.rds'))
  d0_write(list(origin=branch,seed=seed,NP=256L,dimensions=14L,initialization='Prespecified independent 50:50 full-domain/growth-feasibility proposals; no truth/old-fit injection',data_sha256=c$data_sha),file.path(out,'origin.tsv'))
  deadline<-as.numeric(Sys.time())+as.numeric(controls$global_soft_seconds)
  cl<-d0_cluster(s,task,FALSE,deadline,as.integer(controls$global_worker_max_callbacks));on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  native_fn<-function(x).GlobalEnv$d0_ev$fn(x)
  native<-NULL;termination<-'native_return'
  native<-tryCatch(DEoptim::DEoptim(native_fn,rep(0,14),rep(1,14),control=list(NP=256L,
       itermax=as.integer(controls$global_max_generations),trace=25L,cluster=cl,initialpop=population,reltol=1e-7,steptol=300L)),
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
 } else if(startsWith(mode,'local_')||mode=='positive_local') {
  source<-if(positive)file.path(s,'truth/truth_native.rds')else file.path(s,'checks',paste0('globalcheck_',branch),'verified_candidate.rds');initial<-readRDS(source)
  if(positive)initial$data_sha256<-c$data_sha
  stopifnot(identical(initial$data_sha256,c$data_sha));deadline<-as.numeric(Sys.time())+as.numeric(controls$local_soft_seconds)
  ev<-d0_evaluator(c,out,task,deadline,as.integer(controls$local_max_callbacks),infeasible_infinite=TRUE)
  if(!positive)stopifnot(!is.null(initial$unit),identical(d0_parameters(c,initial$unit),initial$params))
  start<-if(positive)d0_unit(c,initial$params)else initial$unit
  q<-ev$fn(start);stopifnot(abs(q-initial$loss$objective)<=1e-8)
  d0_rds(start,file.path(out,'actual_local_start.rds'))
  d0_write(list(pass=TRUE,source_task=if(positive)'separate_truth_control'else paste0('globalcheck_',branch),first_score=q,saved_source_score=initial$loss$objective,difference=q-initial$loss$objective,projection=FALSE,native_unit_preserved=!positive),file.path(out,'start_identity.tsv'))
  sb<-tryCatch(nloptr::nloptr(start,ev$fn,lb=rep(0,14),ub=rep(1,14),opts=list(algorithm='NLOPT_LN_SBPLX',
      maxeval=if(positive)256L else max(2L,as.integer(controls$local_max_callbacks*.7)),maxtime=max(1,(deadline-as.numeric(Sys.time()))*.7),xtol_abs=rep(1e-4,14),ftol_abs=1e-8,print_level=0)),
      d0_budget_stop=function(e)NULL)
  if(!is.null(sb))d0_rds(sb,file.path(out,'native_sbplx.rds'))
  best<-d0_best(out);remaining<-as.integer(controls$local_max_callbacks)-ev$state$n
  pattern<-NULL
  source(file.path(Sys.getenv('D0_RELEASE'),'pattern.R'))
  if(remaining>5L&&as.numeric(Sys.time())<deadline)pattern<-tryCatch(d0_pattern(ev$fn,best$unit,o2review_directions(c$map),remaining,deadline,mesh=c(.02,.005,.001,.0001),tolerance=1e-7),d0_budget_stop=function(e)NULL)
  if(!is.null(pattern))d0_rds(pattern,file.path(out,'native_pattern.rds'))
  winner<-d0_best(out);d0_rds(winner,file.path(out,'accepted_candidate.rds'))
  d0_write(list(pass=TRUE,origin=branch,method='NLOPT_LN_SBPLX_then_bounded_pattern',
       termination=if(ev$state$n>=as.integer(controls$local_max_callbacks))'callback_limit'else if(as.numeric(Sys.time())>=deadline)'time_limit'else'native_return',
       evaluations=ev$state$n,objective=winner$loss$objective,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'))
 } else if(is_check) {
  candidates<-list();proof<-list()
  stages<-if(positive)'positive_local'else if(startsWith(mode,'globalcheck_'))paste0('global_',branch)else c(paste0('global_',branch),paste0('local_',branch))
  for(stage in stages) {
   source<-file.path(s,'fits',stage,'accepted_candidate.rds');if(!file.exists(source))next
   old<-readRDS(source);stopifnot(identical(old$data_sha256,c$data_sha),!is.null(old$unit),identical(d0_parameters(c,old$unit),old$params));runs<-d0_forward(c,old$params,task);pred<-d0_predict(c,runs,old$params);fresh<-d0_loss(c,pred)
   stopifnot(abs(fresh$objective-old$loss$objective)<=1e-8,isTRUE(all.equal(pred,old$predictions,tolerance=1e-12,check.attributes=TRUE)))
   proof[[length(proof)+1L]]<-data.frame(source_task=stage,saved=old$loss$objective,fresh=fresh$objective,difference=fresh$objective-old$loss$objective,independent=fresh$independent)
   candidates[[length(candidates)+1L]]<-list(params=old$params,unit=old$unit,predictions=pred,runs=runs,loss=fresh,source_task=stage,data_sha256=c$data_sha)
  }
  stopifnot(length(candidates)>0L);best<-candidates[[which.min(vapply(candidates,function(z)z$loss$objective,numeric(1)))]]
  d0_rds(best,file.path(out,'verified_candidate.rds'));d0_emit(do.call(rbind,proof),file.path(out,'candidate_checks.tsv'))
  if(startsWith(mode,'globalcheck_')) {
   d0_write(list(pass=TRUE,origin=branch,objective=best$loss$objective,data_sha256=c$data_sha,wall_seconds=proc.time()[3]-started),file.path(out,'result.tsv'));d0_guard(s,out);d0_seal(out);return(invisible(NULL))
  }
  # Fresh directional checks are attached to this exact score-selected vector.
  ev<-d0_evaluator(c,out,task,as.numeric(Sys.time())+3000,60L,save_native=FALSE,infeasible_infinite=TRUE)
  u<-best$unit;stopifnot(!is.null(u),identical(d0_parameters(c,u),best$params));directions<-list()
  for(h in c(.001,.0001))for(i in seq_along(u))for(sign in c(-1,1)) {
   v<-u;v[i]<-min(1,max(0,v[i]+sign*h));if(identical(v,u))next
   q<-ev$fn(v);directions[[length(directions)+1L]]<-data.frame(mesh=h,parameter=c$spec$param_name[i],direction=sign,objective=q,improvement=best$loss$objective-q)
  }
  directions<-do.call(rbind,directions);d0_emit(directions,file.path(out,'directional_checks.tsv'));stable<-all(directions$improvement<=1e-7);directional_calls<-ev$state$n
  # Truth is opened only after score-based estimate selection has finished.
  truth<-readRDS(file.path(s,'truth/truth_parameters.rds'));p<-c$map$table
  p$truth<-vapply(p$parameter,function(n)truth[[n]],numeric(1));p$estimate<-vapply(p$parameter,function(n)best$params[[n]],numeric(1))
  p$natural_error<-p$estimate-p$truth;p$log_ratio<-ifelse(p$truth>0&p$estimate>0,log(p$estimate/p$truth),NA_real_)
  p$scaled_error<-(o2prof_transform(p$estimate,p)-o2prof_transform(p$truth,p))/(o2prof_transform(p$upper,p)-o2prof_transform(p$lower,p));p$estimated<-p$status=='active'
  p$engineering_parameter_resolution<-abs(p$scaled_error)<=.01
  p$boundary_lower<-abs(p$estimate-p$lower)<=64*.Machine$double.eps*pmax(1,abs(p$lower));p$boundary_upper<-abs(p$estimate-p$upper)<=64*.Machine$double.eps*pmax(1,abs(p$upper))
  p$endpoint_data_active<-TRUE
  if(best$params$p_misseg==0)p$endpoint_data_active[p$parameter=='k_o_mis']<-FALSE
  if(best$params$buffer_beta==0)p$endpoint_data_active[p$parameter=='buffer_n_exp']<-FALSE
  if(best$params$buffer_smax==0)p$endpoint_data_active[p$parameter%in%c('buffer_beta','buffer_n_exp')]<-FALSE
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
    estimated_response<-c(lambda=ei$lambda,mu=ei$mu,p_per_chromosome=ei$p_per_chromosome,survival=sv(best$params),effective_error=ee$effective_error_probability_conditional_on_viability,biological_daughter_yield=ee$biological_viable_daughters_per_ordinary_division,retained_loss_rate=ei$ordinary_retained_loss,retained_gain_rate=ei$ordinary_retained_gain)
    response[[length(response)+1L]]<-data.frame(N=N,oxygen_pct=o,function_name=names(tv),truth=unname(tv),estimate=unname(estimated_response),error=unname(estimated_response-tv))
   }
  }
  d0_emit(do.call(rbind,response),file.path(out,'function_recovery.tsv'))
  truth_score<-readRDS(file.path(s,'truth/truth_native.rds'))$loss$objective;delta<-best$loss$objective-truth_score
  stopifnot(delta>=-1e-8)
  gate<-delta<=1e-6&&best$loss$rms_growth_error<=.001&&best$loss$max_growth_error<=.01&&best$loss$max_kary_tv<=.001&&best$loss$max_flow_tv<=.001
  d0_write(c(best$loss,list(pass=TRUE,origin=branch,prediction_recovery=gate,
       local_directional_stability=stable,objective_delta=delta,independent_blind_origin=!positive,
       source_task=best$source_task,model_forward_calls=length(candidates)+directional_calls,candidate_replay_calls=length(candidates),directional_evaluation_calls=directional_calls,data_sha256=c$data_sha,
       all_parameters_unique=FALSE,wall_seconds=proc.time()[3]-started)),file.path(out,'result.tsv'))
  pdf(file.path(out,'parameter_recovery.pdf'),width=11.25,height=7.2)
  active<-p[p$estimated,];par(mar=c(10,5,3,1));barplot(abs(active$scaled_error),names.arg=active$parameter,las=2,col='#3d6da8',ylab='Absolute error / transformed domain width',main=paste('D0 conditional known-nuisance recovery; origin',branch));dev.off()
  pdf(file.path(out,'prediction_recovery.pdf'),width=8.75,height=5.625);barplot(unlist(best$loss[c('max_growth_error','max_kary_tv','max_flow_tv')]),names.arg=c('Growth / true sigma','Karyotype TV','Flow TV'),col='#3d6da8',ylab='Maximum discrepancy',main=paste('Expected observations; origin',branch));dev.off()
 } else stop('Unsupported fit mode')
 d0_guard(s,out);d0_seal(out)
 },error=function(e){d0_write(list(failed=TRUE,message=conditionMessage(e),wall_seconds=proc.time()[3]-started),file.path(out,'FAILED.tsv'));stop(e)})
}
main()
