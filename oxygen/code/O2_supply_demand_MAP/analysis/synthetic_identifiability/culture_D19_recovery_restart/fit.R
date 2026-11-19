# Six bounded population-linked global blocks; D19 only and no truth-based seeding.
d_chain<-c('import','bench','de1b',paste0('de',2:6),'local','pattern','verify1','replay1a','replay1b','verify2','replay2a','replay2b')
d_fit<-function(study,task,out,env,check){
 origin<-substr(task,1,3);stage<-substring(task,5);ix<-match(stage,d_chain);stopifnot(!is.na(ix))
 input<-d_read(study,origin);c<-d_context(env,input$case,input$contract,input$id)
 previous<-if(ix>1L)readRDS(file.path(study,'verification',paste0(origin,'_',d_chain[ix-1]),'state.rds'))else NULL
 seed<-readRDS(file.path(study,'verification/T00/rng',paste0(origin,'.rds')))
 if(!is.null(previous))d_validate_state(c,previous,origin)
 assign('.Random.seed',if(is.null(previous))seed$state else previous$rng,.GlobalEnv)
 ev<-d_evaluator(c,out,origin,stage,previous);on.exit(ev$close(),add=TRUE);complete<-FALSE
 inherited<-if(is.null(previous))list()else previous$quality_records
 finish<-function(extra=list()){
  if(is.null(extra$quality_records))extra$quality_records<-inherited
  r_result(out,ev$snapshot(extra));if(is.finite(ev$state$best$score))sf_emit(data.frame(origin=origin,key=ev$state$best$key,D=ev$state$best$score),file.path(out,'best_identity.tsv'));complete<<-TRUE
 }
 on.exit(if(!complete)r_result(out,ev$snapshot(list(status='FAILED_PARTIAL_ATTEMPTS_CHARGED'))),add=TRUE)
 if(stage=='import'){
  bank<-d_bank(study,origin);u<-do.call(rbind,lapply(bank$members,`[[`,'unit'))
  check('bank_contract_and_shape',bank$contract_sha==sf_hash(c$contract)&&identical(dim(u),c(200L,19L)))
  if(!is.null(c$id$bank_lock))check('locked_native_order',identical(sf_hash(u),c$id$bank_lock$population_sha))
  for(i in 1:200){p<-d_decode(c,as.numeric(u[i,]));check(paste0('member_',i,'_native19'),identical(p[c$free$parameter],bank$members[[i]]$params[c$free$parameter]))}
  sf_rds(list(population=u,rng=.Random.seed),file.path(out,'requested_population.rds'));finish(list(status='READY',population=u));return(invisible(NULL))
 }
 if(stage=='bench'||startsWith(stage,'de')){
  suppressPackageStartupMessages(library(DEoptim));pop<-previous$population;gens<-if(stage=='bench')1L else if(stage=='de1b')59L else 60L
  stopifnot(identical(dim(pop),c(200L,19L)));sf_rds(list(population=pop,rng=.Random.seed),file.path(out,'requested_population.rds'))
  bank<-if(stage=='bench')d_bank(study,origin)else NULL
  fn<-function(u){start<-proc.time()[3];q<-ev$fn(u)
   if(stage=='bench'&&ev$state$n<=200L){
    i<-ev$state$n;check(paste0('imported_physics_',i),is.finite(q)&&ev$state$last$feasible)
    f_compare(ev$state$last$pred,bank$members[[i]]$pred,check,paste0('imported_member_',i))
   };ev$state$times[[ev$state$n]]$seconds<-proc.time()[3]-start;q}
  de<-DEoptim::DEoptim(fn,rep(0,19),rep(1,19),control=list(NP=200L,itermax=gens,trace=FALSE,initialpop=pop,strategy=2L,F=.8,CR=.9,reltol=0,steptol=gens+1L,VTR=-Inf,parallelType=0))
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
   check(paste0('replay_',j,'_native_parameters'),identical(d_decode(c,candidate$unit)[cg_roles],candidate$params))
   v<-d_independent_loss(c,z$pred);check(paste0('replay_',j,'_independent_components'),max(abs(v-c(growth=z$loss$growth,kary=z$loss$kary,value=score)))<=1e-8)
   full<-list(candidate=candidate,pred=z$pred,loss=z$loss,task=task,process=Sys.getpid())
   sf_rds(full,file.path(out,paste0('point_',candidate$key,'.rds')))
   if(which=='b'){
    a<-readRDS(file.path(study,'verification',paste0(origin,'_replay',pass,'a'),paste0('point_',candidate$key,'.rds')))
    check(paste0('two_fresh_processes_',j),a$task!=task)
    f_compare(a$pred,z$pred,check,paste0('two_process_',j));check(paste0('probability_normalization_',j),all(vapply(z$pred$kary,function(p)abs(sum(p)-1)<=1e-12,logical(1))))
    records[[candidate$key]]<-list(candidate=candidate,verified=TRUE,first_task=a$task,second_task=task,state_file=file.path('verification',task,paste0('point_',candidate$key,'.rds')))
   }
  }
  if(which=='b'){
   inherited<-modifyList(inherited,records)
   rows<-lapply(inherited,function(r){p<-r$candidate;local<-previous$local_certificates[[p$key]];data.frame(origin=origin,key=p$key,D=p$score,D_G=p$loss$growth,D_K=p$loss$kary,rms_growth=p$loss$rms_growth,max_growth=p$loss$max_growth,max_kary_tv=p$loss$max_kary_tv,verified=TRUE,prediction_pass=p$score<=1e-4&&p$loss$rms_growth<=.01&&p$loss$max_growth<=.03&&p$loss$max_kary_tv<=.005,checked=!is.null(local),local_pass=if(is.null(local))FALSE else local$complete&&local$max_gain<=1e-5,first_process_task=r$first_task,second_process_task=r$second_task,state_file=r$state_file)})
   sf_emit(do.call(rbind,rows),file.path(out,'quality.tsv'))
  }
  finish(list(status=if(stage=='replay2b')'COMPLETE'else'READY',replay_pending=pending,local_certificates=previous$local_certificates,quality_records=inherited));return(invisible(NULL))
 }
 x<-previous$best$unit;stopifnot(length(x)==19,is.finite(previous$best$score))
 replay<-ev$fn(x,bucket='handoff');check('exact_native_handoff_score',abs(replay-previous$best$score)<=1e-8);check('exact_native_handoff_parameters',identical(d_decode(c,x)[cg_roles],previous$best$params))
 if(stage=='local'){
  suppressPackageStartupMessages(library(nloptr));stopped<-FALSE
  nl<-tryCatch(nloptr::nloptr(x,function(u)ev$fn(u),lb=rep(0,19),ub=rep(1,19),opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=6000L,maxtime=0,xtol_abs=rep(1e-5,19),ftol_abs=1e-8,print_level=0)),recovery_budget_stop=function(e){stopped<<-TRUE;NULL})
  sf_rds(nl,file.path(out,'native_Subplex_result.rds'));finish(list(status='READY',wrapper_allowance_exhausted=stopped,native_status=if(is.null(nl))NA_integer_ else nl$status));return(invisible(NULL))
 }
 if(stage=='pattern'){
  pol<-r_poll(ev$fn,x,r_directions(c),c(.02,.005,.001,.0001),1e-8,2000L-ev$state$counts['pattern']);sf_emit(pol$rows,file.path(out,'pattern_polls.tsv'));finish(list(status='READY',pattern_budget_limited=pol$budget_limited));return(invisible(NULL))
 }
 stopifnot(stage%in%c('verify1','verify2'))
 anchor<-d_candidate(c,x,ev$state$last,origin,stage);directions<-r_directions(c);rows<-list();seen<-character();all_done<-TRUE;gain<-0
 for(h in c(.001,.0001,.00002))for(j in seq_len(nrow(directions)))for(sign in c(-1,1)){
  u<-r_admissible(x,sign*directions[j,],h);if(is.null(u))next;k<-sf_hash(u);if(k%in%seen)next;seen<-c(seen,k)
  q<-tryCatch(ev$fn(u),recovery_budget_stop=function(e)NA_real_);if(is.na(q)){all_done<-FALSE;next};gain<-max(gain,replay-q)
  rows[[length(rows)+1L]]<-data.frame(anchor_key=anchor$key,mesh=h,direction=j,sign=sign,actual_step=sqrt(sum((u-x)^2)),key=k,score=q,improvement=replay-q)
 }
 sf_emit(if(length(rows))do.call(rbind,rows)else data.frame(anchor_key=character()),file.path(out,'verification_polls.tsv'))
 cert<-if(is.null(previous$local_certificates))list()else previous$local_certificates;cert[[anchor$key]]<-list(complete=all_done,max_gain=gain,task=task)
 pending<-list(anchor,ev$state$best);if(!is.null(ev$state$gate))pending<-c(pending,list(ev$state$gate));pending<-pending[!duplicated(vapply(pending,`[[`,character(1),'key'))]
 sf_rds(anchor,file.path(out,'checked_representative.rds'));finish(list(status='READY',checked_representative=anchor,replay_pending=pending,local_certificates=cert))
}
