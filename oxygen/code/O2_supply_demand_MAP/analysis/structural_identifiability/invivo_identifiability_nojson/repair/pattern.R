# Bounded normalized block/pattern search with independently checked observation moves.
pattern_directions <- function(names_active) {
 blocks<-list(resource=c('log10_o2_S0','log10_kappa_O','log10_eta_o2'),growth=c('log10_lam_max','log10_alpha_o2','gamma_growth'),
  death=c('log10_mu_hp','gamma_mu','log10_O2_crit','n_O'),chromosome=c('log10_p_mis_base','log10_p_misseg','log10_k_o_mis','log10_p_wgd','buffer_smax','log10_buffer_beta','log10_buffer_n_exp'),clearance='log10_k_clear')
 dirs<-list()
 for(n in names_active) for(sign in c(-1,1)) {v<-setNames(rep(0,length(names_active)),names_active);v[n]<-sign;dirs[[paste0(n,'_',sign)]]<-v}
 for(block in names(blocks)) {
  active<-intersect(blocks[[block]],names_active);if(length(active)<2L) next
  for(kind in c('same','alternating')) for(sign in c(-1,1)) {
   v<-setNames(rep(0,length(names_active)),names_active);v[active]<-sign*if(kind=='same') rep(1,length(active)) else rep(c(1,-1),length.out=length(active))
   dirs[[paste(block,kind,sign,sep='_')]]<-v
  }
 }
 dirs
}
run_pattern_local <- function(study,task,spec,fixture=FALSE) {
 out<-file.path(study,if(fixture) 'tests' else 'fits',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 stopifnot(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))==1L)
 ctx<-formal_context(study,spec);fixed<-spec$fixed;start<-resolve_start(ctx,spec);p<-start$parameters
 if(length(fixed))p[names(fixed)]<-fixed
 atomic_rds(spec,file.path(out,'actual_spec.rds'));atomic_rds(list(parameters=p,resolution=start),file.path(out,'actual_local_start.rds'))
 ctx$provenance$start_sha256<-hash_file(file.path(out,'actual_local_start.rds'))
 free<-setdiff(names(p),names(fixed));dynamic<-setdiff(free,c('log10_rho_2N','log10_sigma_burden'));span<-ctx$upper-ctx$lower
 .ivi_evaluator<<-new_evaluator(ctx,out,fixed,'bounded_block_pattern',getOption('ivi.job.deadline'))
 t<-proc.time()[[3]];start_value<-.ivi_evaluator(p[free])
 if(is.null(fixed) && !is.null(start$expected))stopifnot(abs(start_value-start$expected)<1e-8)
 events<-list();proposals<-list();completed_scales<-numeric();limited<-FALSE;call_limited<-FALSE
 maxcalls<-if(fixture) 160L else 10000L;maxsweeps<-if(fixture) 1L else 25L
 evaluate<-function(q,kind,step=NA_real_,direction='') {
  if(.ivi_state$n>=maxcalls)stop('IVI_CALL_LIMIT')
  value<-.ivi_evaluator(q[free]);events[[length(events)+1L]]<<-list(parameters=q,objective=value,components=.ivi_state$last_components,kind=kind,step=step,direction=direction,evaluation=.ivi_state$n)
  value
 }
 observe<-function() {
  prior<-.ivi_state$best_record;proposal<-observation_proposal(ctx,prior,fixed)
  if(max(abs(proposal$parameters-prior$parameters))<1e-14)return(invisible(NULL))
  val<-evaluate(proposal$parameters,'observation');gap<-val-proposal$expected
  proposals[[length(proposals)+1L]]<<-c(proposal,list(actual=val,difference=gap))
  if(abs(gap)>1e-8)stop('Observation-block canonical equivalence failure')
  invisible(NULL)
 }
 stop_reason<-'complete'
 tryCatch({
  observe()
  directions<-pattern_directions(dynamic)
  for(step in if(fixture) c(.001) else c(.02,.005,.001,.0001)) {
   for(sweep in seq_len(maxsweeps)) {
    base<-.ivi_state$best_record$parameters;before<-.ivi_state$best
    for(name in names(directions)) {
     q<-base;q[dynamic]<-base[dynamic]+step*span[dynamic]*directions[[name]]
     if(any(q<ctx$lower | q>ctx$upper))next
     evaluate(q,'poll',step,name)
    }
    if(.ivi_state$best<before-1e-10)observe()
    if(.ivi_state$best>=before-1e-10)break
   }
   completed_scales<-c(completed_scales,step)
  }
 },error=function(e) {
  msg<-conditionMessage(e)
  if(grepl('IVI_SOFT_TIME_LIMIT',msg,fixed=TRUE)){limited<<-TRUE;stop_reason<<-'soft_time_limit'}
  else if(grepl('IVI_CALL_LIMIT',msg,fixed=TRUE)){call_limited<<-TRUE;stop_reason<<-'call_limit'}
  else stop(e)
 })
 best<-incumbents(out,ctx);checked<-validate_point(ctx,best$parameters,record=best)
 stopifnot(checked$components$L<=start_value+1e-8,all(best$parameters[names(fixed)]==fixed))
 atomic_rds(events,file.path(out,'polish_evaluations.rds'));atomic_rds(proposals,file.path(out,'observation_proposals.rds'))
 atomic_rds(list(method='bounded_block_pattern',termination=stop_reason,completed_scales=completed_scales,max_sweeps_per_scale=maxsweeps,max_calls=maxcalls,par=best$parameters,value=best$components$L),file.path(out,'native_result.rds'))
 atomic_rds(list(record=best,checked=checked$components,predictions_full=checked$predictions,fixed=fixed,policy=spec$policy,budget_limited=limited,call_limited=call_limited,fixture=fixture,provenance=ctx$provenance),file.path(out,'accepted_candidate.rds'))
 for(k in names(checked$predictions$tables))write_tsv(checked$predictions$tables[[k]],file.path(out,paste0(k,'_predictions.tsv')))
 write_kv(list(pass=TRUE,policy=spec$policy,method='bounded_block_pattern',termination=stop_reason,start_objective=start_value,objective=checked$components$L,wall_seconds=proc.time()[[3]]-t,
  evaluations=.ivi_state$n,objective_seconds=.ivi_state$seconds,improvements=.ivi_state$improvements,budget_limited=limited,call_limited=call_limited,
  native_convergence=NA_integer_,fixed_residual=if(length(fixed))max(abs(best$parameters[names(fixed)]-fixed)) else 0,
  observation_proposals=length(proposals),free_coordinates=length(free),dynamic_coordinates=length(dynamic),workers=1L),file.path(out,'result.tsv'))
 write_kv(meter_values(),file.path(out,'runtime_meter.tsv'));finalize_payload(out)
}
