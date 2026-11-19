# Finite feasibility-only construction. No D19 or noisy response is read/scored.
f_write_bank<-function(c,out,origin,state){
 stopifnot(length(state$archive)>=200L,length(state$anchors)==8L)
 idx<-c(state$anchors,head(setdiff(seq_along(state$archive),state$anchors),192L));stopifnot(length(idx)==200L)
 members<-state$archive[idx];stopifnot(all(vapply(members,function(a)!is.null(a$pred$states),logical(1))))
 bank<-list(origin=origin,members=members,contract_sha=sf_hash(c$contract),rng=state$rng,construction_attempts=state$count,kind='truth_independent_conditional_feasible_initial_population_not_recovery')
 sf_rds(bank,file.path(out,'bank.rds'));u<-do.call(rbind,lapply(members,`[[`,'unit'));colnames(u)<-f_free(c$d)$parameter
 sf_emit(data.frame(row=1:200,key=vapply(members,`[[`,character(1),'key'),is_anchor=1:200<=8,source_call=vapply(members,`[[`,integer(1),'call'),u),file.path(out,'bank_members.tsv'))
 sf_emit(data.frame(parameter=colnames(u),unit_min=apply(u,2,min),unit_max=apply(u,2,max),span=apply(u,2,function(v)diff(range(v))),lower_hits=colSums(u==0),upper_hits=colSums(u==1)),file.path(out,'bank_diversity.tsv'))
}
f_build<-function(study,task,out,env,check){
 origin<-sub('_[^_]+$','',task);stage<-sub('^.*_','',task);chain<-c('bench','de1b',paste0('de',2:6),'fill');position<-match(stage,chain);stopifnot(!is.na(position))
 contract<-readRDS(file.path(study,'inputs/public/contract.rds'));c<-f_context(env,contract);rng<-readRDS(file.path(study,'verification/T00/rng',paste0(origin,'.rds')))
 previous<-if(position==1)NULL else readRDS(file.path(study,'verification',paste0(origin,'_',chain[position-1]),'state.rds'))
 assign('.Random.seed',if(is.null(previous))rng$state else previous$rng,.GlobalEnv)
 ev<-f_evaluator(c,out,'constructor',origin,previous);on.exit(ev$close(),add=TRUE);complete<-FALSE
 on.exit(if(!complete)f_save(out,ev$snapshot(list(status='FAILED_PARTIAL'))),add=TRUE)
 target<-FALSE;cap<-FALSE;pop<-if(is.null(previous))matrix(runif(3800),200,19)else previous$population;native<-NULL
 if(stage!='fill'){
  generations<-if(stage=='bench')2L else if(stage=='de1b')6L else 8L;stopifnot(identical(dim(pop),c(200L,19L)))
  sf_rds(list(population=pop,rng=.Random.seed),file.path(out,'requested_population.rds'))
  suppressPackageStartupMessages(library(DEoptim));native<-tryCatch(DEoptim::DEoptim(ev$fn,rep(0,19),rep(1,19),control=list(NP=200L,itermax=generations,initialpop=pop,trace=FALSE,strategy=2L,F=.8,CR=.9,reltol=0,steptol=generations+1L,VTR=-Inf,parallelType=0)),feasibility_bank_complete=function(e){target<<-TRUE;NULL},feasibility_budget_stop=function(e){cap<<-TRUE;NULL})
  actual<-do.call(rbind,ev$state$initial);sf_rds(actual,file.path(out,'actual_initial_population.rds'))
  check('actual_native_population_prefix',identical(unname(actual),unname(pop[seq_len(nrow(actual)),,drop=FALSE])))
  if(!is.null(native)){check('native_full_callback_count',native$optim$nfeval==ev$state$n&&ev$state$n==200L*(generations+1L));pop<-native$member$pop}
  sf_rds(native,file.path(out,'native_DE_result.rds'))
  sf_emit(data.frame(call=seq_len(ev$state$n),component=ifelse(seq_len(ev$state$n)<=200,if(stage=='bench')'full_box_uniform_initial'else'linked_initial_population','native_DE_strategy2_trial')),file.path(out,'proposal_components.tsv'))
 }else{
  components<-list();initial_count<-ev$state$counts
  while(ev$state$counts<12000L){
   call<-ev$state$counts-initial_count+1L;radii<-c(.25,.10,.03);h<-radii[(call-1L)%%3L+1L];k<-NA_integer_
   if(runif(1)<.25||!length(ev$state$anchors)){u<-runif(19);component<-'full_box_uniform'}else{
    k<-ev$state$anchors[(call-1L)%%length(ev$state$anchors)+1L];u<-ev$state$archive[[k]]$unit+runif(19,-h,h);u<-1-abs((u%%2)-1);component<-'own_feasible_anchor_reflected_perturbation'
   }
   components[[call]]<-data.frame(call=ev$state$counts+1L,component=component,radius=if(is.na(k))NA_real_ else h,anchor_archive_index=k)
   tryCatch(ev$fn(u),feasibility_bank_complete=function(e){target<<-TRUE},feasibility_budget_stop=function(e){cap<<-TRUE});if(target||cap)break
  }
  sf_emit(do.call(rbind,components),file.path(out,'proposal_components.tsv'))
 }
 state<-ev$snapshot(list(status=if(target)'BANK_READY'else if(stage=='fill'||cap)'FEASIBILITY_CONSTRUCTION_UNRESOLVED'else'CONTINUE',population=pop,native_completed=!is.null(native),bank_target_reached=target,partial_native_return=target||cap))
 if(target)f_write_bank(c,out,origin,state)
 f_save(out,state);complete<-TRUE
}
