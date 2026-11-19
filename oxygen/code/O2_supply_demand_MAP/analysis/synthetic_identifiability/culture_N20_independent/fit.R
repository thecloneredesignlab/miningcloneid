# Six bounded population-linked global blocks; N20 only and no truth-based seeding.
n_chain<-c('import','bench','de1b',paste0('de',2:6),'local','pattern','verify1','replay1a','replay1b','verify2','replay2a','replay2b')
n_fit<-function(study,task,out,env,check,phase_name=task){
 origin<-substr(task,1,3);stage<-substring(task,5);ix<-match(stage,n_chain);stopifnot(!is.na(ix))
 input<-n_data(study,origin);c<-n_context(env,input$case,input$contract,input$id)
 previous<-if(ix>1L)readRDS(file.path(study,'verification',paste0(origin,'_',n_chain[ix-1]),'state.rds'))else NULL
 if(!is.null(previous))stopifnot(identical(previous$stage,n_chain[ix-1]))
 seed<-readRDS(file.path(study,'verification/T00/rng',paste0(origin,'.rds')))
 stopifnot(identical(sf_hash(seed$state),input$id$rng_initial_sha))
 if(!is.null(previous))n_validate_state(c,previous,origin)
 assign('.Random.seed',if(is.null(previous))seed$state else previous$rng,.GlobalEnv)
 ev<-n_evaluator(c,out,origin,stage,previous);on.exit(ev$close(),add=TRUE);complete<-FALSE
 inherited<-if(is.null(previous))list()else previous$quality_records
 finish<-function(extra=list()){
  if(is.null(extra$quality_records))extra$quality_records<-inherited
  r_result(out,ev$snapshot(extra));if(is.finite(ev$state$best$score))sf_emit(data.frame(origin=origin,key=ev$state$best$key,Q=ev$state$best$score),file.path(out,'best_identity.tsv'));complete<<-TRUE
 }
 on.exit(if(!complete)r_result(out,ev$snapshot(list(status='FAILED_PARTIAL_ATTEMPTS_CHARGED'))),add=TRUE)
 if(stage=='import'){
  bank<-n_bank(study,origin);u<-bank$population
  check('bank_shape_order_and_roles',identical(dim(u),c(200L,20L)) &&
        identical(bank$roles,c$free$parameter) && bank$bounds_sha==sf_hash(c$free) &&
        identical(sf_hash(u),c$id$population_sha))
  check('no_responses_or_legacy_checkpoints',identical(names(bank),
        c('origin','source_sha','bounds_sha','roles','population','physical_native')))
  for(i in 1:200){
   p<-n_decode(c,as.numeric(u[i,]))
   check(paste0('member_',i,'_native20'),
         identical(p[setdiff(c$free$parameter,'sigma_growth')],bank$physical_native[[i]]) &&
         p$sigma_growth>0)
  }
  sf_rds(list(population=u,rng=.Random.seed),file.path(out,'requested_population.rds'));finish(list(status='READY',population=u));return(invisible(NULL))
 }
 if(stage=='bench'||startsWith(stage,'de')){
  suppressPackageStartupMessages(library(DEoptim));pop<-previous$population;gens<-if(stage=='bench')1L else if(stage=='de1b')59L else 60L
  stopifnot(identical(dim(pop),c(200L,20L)));sf_rds(list(population=pop,rng=.Random.seed),file.path(out,'requested_population.rds'))
  fn<-function(u){start<-proc.time()[3];q<-ev$fn(u)
   if(stage=='bench'&&ev$state$n<=200L){
    i<-ev$state$n;check(paste0('imported_physics_',i),is.finite(q)&&ev$state$last$feasible)
   };ev$state$times[[ev$state$n]]$seconds<-proc.time()[3]-start;q}
  check('native_block_fits_remaining_attempts',ev$state$counts['global']+200L*(gens+1L)<=n_limits['global'])
  de<-DEoptim::DEoptim(fn,rep(0,20),rep(1,20),control=list(NP=200L,itermax=gens,trace=FALSE,initialpop=pop,strategy=2L,F=.8,CR=.9,reltol=0,steptol=gens+1L,VTR=-Inf,parallelType=0))
  actual<-do.call(rbind,ev$state$initial);check('exact_native_initial_order',identical(unname(actual),unname(pop)));check('native_callback_count',de$optim$nfeval==ev$state$n&&ev$state$n==200L*(gens+1L))
  sf_rds(actual,file.path(out,'actual_initial_population.rds'));sf_rds(de,file.path(out,'native_DE_result.rds'))
  finish(list(status='READY',population=de$member$pop,native_evaluations=de$optim$nfeval,generations=gens));return(invisible(NULL))
 }
 if(startsWith(stage,'replay')){
  pass<-as.integer(substr(stage,7,7));which<-substr(stage,8,8)
  # Both a and b are separate allocations/R processes, with unchanged exact vectors.
  pending<-previous$replay_pending;stopifnot(length(pending)>=1)
  records<-list()
  for(j in seq_along(pending)){
   candidate<-pending[[j]];score<-ev$fn(candidate$unit);z<-ev$state$last
   check(paste0('replay_',j,'_score'),z$feasible&&abs(score-candidate$score)<=1e-8)
   check(paste0('replay_',j,'_native_parameters'),identical(n_decode(c,candidate$unit)[cg_roles],candidate$params))
   v<-n_independent_loss(c,z$pred,n_decode(c,candidate$unit));check(paste0('replay_',j,'_independent_components'),max(abs(v-c(growth=z$loss$growth,kary=z$loss$kary,value=score)))<=1e-8)
   full<-list(candidate=candidate,pred=z$pred,loss=z$loss,task=phase_name,process=Sys.getpid(),slurm_job=Sys.getenv('SLURM_JOB_ID'))
   sf_rds(full,file.path(out,paste0('point_',candidate$key,'.rds')))
   if(which=='b'){
    a<-readRDS(file.path(study,'verification',paste0(origin,'_replay',pass,'a'),paste0('point_',candidate$key,'.rds')))
    check(paste0('two_fresh_processes_',j),a$task!=phase_name && nzchar(full$slurm_job) && a$slurm_job!=full$slurm_job)
    f_compare(a$pred,z$pred,check,paste0('two_process_',j));check(paste0('probability_normalization_',j),all(vapply(z$pred$kary,function(p)abs(sum(p)-1)<=1e-12,logical(1))))
    check(paste0('endpoint_protocol_identity_',j),identical(candidate$pred$selection,z$pred$selection))
    check(paste0('replayed_Q_components_',j),
      max(abs(unlist(a$loss[c('growth','kary','value')])-unlist(z$loss[c('growth','kary','value')])))<=1e-8)
    records[[candidate$key]]<-list(candidate=candidate,verified=TRUE,first_task=a$task,second_task=phase_name,state_file=file.path('verification',task,paste0('point_',candidate$key,'.rds')))
   }
  }
  if(which=='b'){
   inherited<-modifyList(inherited,records)
   rows<-lapply(inherited,function(r){
    p<-r$candidate
    local<-previous$local_certificates[[p$key]]
    data.frame(origin=origin,key=p$key,Q=p$score,growth_log_mean=p$loss$growth,karyotype_log_mean=p$loss$kary,
      verified=TRUE,checked=!is.null(local),local_pass=if(is.null(local))FALSE else
      local$complete&&local$max_gain<=1e-5&&local$variance_gain<=1e-5,
      maximum_gain=if(is.null(local))NA_real_ else local$max_gain,
      variance_gain=if(is.null(local))NA_real_ else local$variance_gain,
      first_process_task=r$first_task,second_process_task=r$second_task,state_file=r$state_file)
   })
   sf_emit(do.call(rbind,rows),file.path(out,'quality.tsv'))
  }
  finish(list(status=if(stage=='replay2b')'COMPLETE'else'READY',replay_pending=pending,local_certificates=previous$local_certificates,quality_records=inherited));return(invisible(NULL))
 }
 x<-previous$best$unit;stopifnot(length(x)==20,is.finite(previous$best$score))
 replay<-ev$fn(x,bucket='handoff');check('exact_native_handoff_score',abs(replay-previous$best$score)<=1e-8);check('exact_native_handoff_parameters',identical(n_decode(c,x)[cg_roles],previous$best$params))
 if(stage=='local'){
  suppressPackageStartupMessages(library(nloptr));stopped<-FALSE
  nl<-tryCatch(nloptr::nloptr(x,function(u)ev$fn(u),lb=rep(0,20),ub=rep(1,20),opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=6000L,maxtime=0,xtol_abs=rep(1e-5,20),ftol_abs=1e-6,print_level=0)),recovery_budget_stop=function(e){stopped<<-TRUE;NULL})
  sf_rds(nl,file.path(out,'native_Subplex_result.rds'));finish(list(status='READY',wrapper_allowance_exhausted=stopped,native_status=if(is.null(nl))NA_integer_ else nl$status));return(invisible(NULL))
 }
 if(stage=='pattern'){
  pol<-r_poll(ev$fn,x,r_directions(c),c(.02,.005,.001,.0001),1e-8,2000L-ev$state$counts['pattern']);sf_emit(pol$rows,file.path(out,'pattern_polls.tsv'));finish(list(status='READY',pattern_budget_limited=pol$budget_limited));return(invisible(NULL))
 }
 stopifnot(stage%in%c('verify1','verify2'))
 anchor<-n_candidate(c,x,ev$state$last,origin,stage);directions<-r_directions(c);rows<-list();seen<-character();all_done<-TRUE;gain<-0
 for(h in c(.001,.0001,.00002))for(j in seq_len(nrow(directions)))for(sign in c(-1,1)){
  u<-r_admissible(x,sign*directions[j,],h);if(is.null(u))next;k<-sf_hash(u);if(k%in%seen)next;seen<-c(seen,k)
  q<-tryCatch(ev$fn(u),recovery_budget_stop=function(e)NA_real_);if(is.na(q)){all_done<-FALSE;next};gain<-max(gain,replay-q)
  rows[[length(rows)+1L]]<-data.frame(anchor_key=anchor$key,mesh=h,direction=j,sign=sign,actual_step=sqrt(sum((u-x)^2)),key=k,score=q,improvement=replay-q)
 }
 sf_emit(if(length(rows))do.call(rbind,rows)else data.frame(anchor_key=character()),file.path(out,'verification_polls.tsv'))
 sigma<-n_sigma_check(c,anchor)
 if(!identical(sigma$unit,anchor$unit)){
  actual<-ev$fn(sigma$unit)
  check('actual_variance_stationarity_score',is.finite(actual)&&abs(actual-sigma$score)<=1e-8)
  check('actual_variance_stationarity_same_means',
     identical(ev$state$last$pred$selection,anchor$pred$selection) &&
     max(abs(ev$state$last$pred$growth-anchor$pred$growth))<=1e-10 &&
     abs(ev$state$last$loss$kary-anchor$loss$kary)<=1e-10)
 }else actual<-anchor$score
 variance_gain<-max(0,anchor$score-actual)
 sf_emit(data.frame(anchor_key=anchor$key,raw_residual_RMS=sigma$unclipped,check_sigma=sigma$scale,
    analytic_Q=sigma$score,replayed_Q=actual,improvement=variance_gain),file.path(out,'variance_stationarity.tsv'))
 cert<-if(is.null(previous$local_certificates))list()else previous$local_certificates
 cert[[anchor$key]]<-list(complete=all_done,max_gain=gain,variance_gain=variance_gain,task=task)
 pending<-list(anchor,ev$state$best)
 pending<-pending[!duplicated(vapply(pending,`[[`,character(1),'key'))]
 sf_rds(anchor,file.path(out,'checked_representative.rds'));finish(list(status='READY',checked_representative=anchor,replay_pending=pending,local_certificates=cert))
}
