# Reproduce each old vector and distinguish both cohort histories, without scoring.
f_diagnose<-function(study,task,out,env,check){
 contract<-readRDS(file.path(study,'verification/P0/feasibility_contract.rds'));c<-f_context(env,contract)
 units<-readRDS(file.path(study,'verification/P0/old_units.rds'));chunk<-as.integer(sub('P2_','',task));ids<-((chunk-1L)*250L+1L):(chunk*250L)
 reference_ids<-unique(as.integer(round(seq(1,1000,length.out=64))));refmode<-FALSE;capture<-FALSE
 backend<-function(u){p<-f_decode(env,c$d,u);if(refmode){z<-f_reference(c,p);z$F<-if(z$feasible)0 else 1.1;z}else f_physics(c,p,capture)}
 ev<-f_evaluator(c,out,'diagnostic',task,backend=backend);on.exit(ev$close(),add=TRUE);complete<-FALSE
 on.exit(if(!complete)f_save(out,ev$snapshot(list(status='FAILED_PARTIAL'))),add=TRUE)
 summary<-list();ledger<-list();references<-list();daily<-list()
 for(i in ids){
  refmode<-FALSE;capture<-i%in%reference_ids;ev$fn(units[i,]);z<-ev$state$last
  check(paste0('old_infeasible_',i),!z$feasible)
  q<-f_short(z);q$candidate<-i;summary[[length(summary)+1L]]<-q;ledger[[length(ledger)+1L]]<-f_ledger(z,i,contract)
  for(cohort in names(z$lineages))for(j in seq_along(z$lineages[[cohort]]$edges)){
   e<-z$lineages[[cohort]]$edges[[j]];daily[[length(daily)+1L]]<-data.frame(candidate=i,cohort=cohort,edge=j,index=seq_along(e$live),day=e$days,reported_live=e$live,state_total=e$row_totals)
  }
  failures<-z$lineages[!vapply(z$lineages,`[[`,logical(1),'feasible')];first<-failures[[1]]$error
  if(capture){
   sf_rds(list(unit=as.numeric(units[i,]),failures=lapply(failures,function(a)list(condition=a$error,parent_state=tail(a$edges,1)[[1]]))),file.path(out,'selected_states',sprintf('%04d.rds',i)))
   refmode<-TRUE;ev$fn(units[i,]);r<-ev$state$last;check(paste0('original_first_failure_',i),!r$feasible&&identical(r$error,first))
   references[[length(references)+1L]]<-data.frame(candidate=i,status='PASS',cohort=first$cohort,parent=first$parent,child=first$child)
  }
 }
 sf_emit(do.call(rbind,summary),file.path(out,'candidate_failure_ledger.tsv'));sf_emit(do.call(rbind,ledger),file.path(out,'passage_supply_ledger.tsv'));sf_emit(do.call(rbind,daily),file.path(out,'daily_supply.tsv'));sf_emit(do.call(rbind,references),file.path(out,'reference_replay_checks.tsv'))
 f_save(out,ev$snapshot(list(status='DIAGNOSTIC_PASS',candidates=length(ids),references=length(references))));complete<-TRUE
}
