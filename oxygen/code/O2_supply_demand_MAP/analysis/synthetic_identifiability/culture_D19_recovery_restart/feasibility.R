# Observation-free physical feasibility and durable candidate accounting.
f_limits<-c(tests=128L,diagnostic=1064L,constructor=12000L,validation=2752L)
f_fail<-function(message)stop(message,call.=FALSE)
f_free<-function(d)d$bounds[d$bounds$parameter!='sigma_growth',,drop=FALSE]
f_decode<-function(env,d,u){
 t<-f_free(d);stopifnot(is.numeric(u),length(u)==19,all(is.finite(u)),all(u>=0&u<=1))
 neutral<-cg_neutral(env,d$bounds);low<-env$o2prof_transform(t$lower,t);high<-env$o2prof_transform(t$upper,t)
 z<-setNames(vapply(seq_len(nrow(d$bounds)),function(i)env$o2prof_transform(neutral[[d$bounds$parameter[i]]],d$bounds[i,,drop=FALSE]),numeric(1)),d$bounds$coordinate)
 z[t$coordinate]<-low+u*(high-low);cg_from_free(env,d,z)
}
f_reduced_condition<-function(e){
 list(class=class(e),message=conditionMessage(e),parent=e$parent_segment$segment_id,child=e$child_segment$segment_id,cohort=e$child_segment$cohort,selection=e$selection,ordinal=e$segment_ordinal,count=e$segment_count)
}
f_sim_valid<-function(sim){
 x<-sim$live_state_obs;stopifnot(is.matrix(x),ncol(x)==133,all(is.finite(x)),min(x)>=0,all(is.finite(sim$Ntot_live_obs)))
}
f_context<-function(env,contract){
 c<-new.env(parent=emptyenv());c$env<-env;c$contract<-contract;c$d<-contract$design;c$hash<-sf_hash(contract);c$core<-env$build_model_core(cfg=c$d$cfg)
 c$extract<-env$ivt_extract_passage_end_state;c$select<-env$ivt_select_segment_observation;c$current<-NULL
 env$ivt_select_segment_observation<-function(...){r<-c$select(...);if(!is.null(c$current)){if(isTRUE(c$current$inside))c$current$target<-r else c$current$scored[[length(c$current$scored)+1L]]<-r};r}
 env$ivt_extract_passage_end_state<-function(sim,reseed_live_cells,grid_pre,target_live_cells=NA_real_,obs_days_local=NULL){
  stopifnot(!is.null(c$current));f_sim_valid(sim);c$current$inside<-TRUE;on.exit(c$current$inside<-FALSE)
  r<-c$extract(sim,reseed_live_cells,grid_pre,target_live_cells,obs_days_local)
  a<-list(selection=r,target=c$current$target,live=as.numeric(sim$Ntot_live_obs),row_totals=rowSums(sim$live_state_obs),days=obs_days_local)
  if(isTRUE(c$current$capture))a$states<-sim$live_state_obs
  c$current$edges[[length(c$current$edges)+1L]]<-a;r
 }
 c
}
f_one<-function(c,a,p,capture=FALSE){
 obs<-new.env(parent=emptyenv());obs$inside<-FALSE;obs$target<-NULL;obs$edges<-list();obs$scored<-list();obs$capture<-capture;c$current<-obs
 on.exit(c$current<-NULL)
 error<-NULL;run<-tryCatch(c$env$ivt_run_lineage(a,cfg=c$d$cfg,run_params=p,model_core=c$core),invitro_protocol_infeasible=function(e){error<<-f_reduced_condition(e);NULL})
 expected<-c$contract$edges[c$contract$edges$cohort==a$cohort & c$contract$edges$kind=='explicit_parent',,drop=FALSE]
 if(is.null(error)){
  stopifnot(!is.null(run),length(obs$edges)==nrow(expected));for(s in run$segment_results)f_sim_valid(s$sim)
  return(list(feasible=TRUE,v=0,passed=nrow(expected),required=nrow(expected),edges=obs$edges,scored=obs$scored,run=run,error=NULL))
 }
 n<-length(obs$edges);stopifnot(n>=1,n<=nrow(expected),error$cohort==a$cohort,error$parent==expected$parent[n],error$child==expected$child[n])
 s<-error$selection;req<-as.numeric(s$required_cells);available<-as.numeric(s$available_cells)
 stopifnot(length(req)==1,is.finite(req),req>0,length(available)==1,is.finite(available),available>=0,available<req)
 stopifnot(!isTRUE(s$passage_executed),sum(vapply(obs$edges,function(e)isTRUE(e$selection$passage_executed),logical(1)))==n-1L)
 r<-available/req;v<-(2*(nrow(expected)-(n-1L))+(1-r))/(2*nrow(expected)+1)
 list(feasible=FALSE,v=v,passed=n-1L,required=nrow(expected),edges=obs$edges,scored=obs$scored,run=NULL,error=error)
}
f_physics<-function(c,p,capture=FALSE){
 z<-lapply(c$d$adapters,function(a)f_one(c,a,p,capture));names(z)<-names(c$d$adapters)
 stopifnot(setequal(names(z),c('2N','4N')));v<-vapply(z,`[[`,numeric(1),'v');ok<-all(v==0)
 pred<-NULL;if(ok){pred<-cg_predict(c$env,c$d,lapply(z,`[[`,'run'),p);stopifnot(length(pred$growth)==114,all(is.finite(pred$growth)),length(pred$kary)==12)}
 list(feasible=ok,F=max(v)+.1*mean(v),lineages=z,pred=pred)
}
f_reference<-function(c,p){
 env<-c$env;wrapped_extract<-env$ivt_extract_passage_end_state;wrapped_select<-env$ivt_select_segment_observation
 env$ivt_extract_passage_end_state<-c$extract;env$ivt_select_segment_observation<-c$select
 on.exit({env$ivt_extract_passage_end_state<-wrapped_extract;env$ivt_select_segment_observation<-wrapped_select})
 error<-NULL;r<-tryCatch(lapply(c$d$adapters,function(a)env$ivt_run_lineage(a,cfg=c$d$cfg,run_params=p,model_core=c$core)),invitro_protocol_infeasible=function(e){error<<-f_reduced_condition(e);NULL})
 list(feasible=is.null(error),error=error,pred=if(is.null(error))cg_predict(env,c$d,r,p)else NULL)
}
f_compare<-function(a,b,check,label){
 check(paste0(label,'_decisions'),identical(a$selection,b$selection))
 check(paste0(label,'_state_identity'),identical(names(a$states),names(b$states)))
 check(paste0(label,'_growth'),identical(names(a$growth),names(b$growth))&&max(abs(a$growth-b$growth))<=1e-10)
 check(paste0(label,'_states'),all(vapply(names(a$states),function(n){x<-a$states[[n]];y<-b$states[[n]];identical(dim(x),dim(y))&&all(abs(x-y)<=1e-10*pmax(1,abs(x),abs(y)))},logical(1))))
 check(paste0(label,'_probability'),max(vapply(names(a$kary),function(n)sum(abs(a$kary[[n]]-b$kary[[n]]))/2,numeric(1)))<=1e-10)
}
f_short<-function(z){
 do.call(rbind,lapply(names(z$lineages),function(cohort){a<-z$lineages[[cohort]];s<-if(is.null(a$error))NULL else a$error$selection
  data.frame(cohort=cohort,feasible=a$feasible,v=a$v,passed=a$passed,required=a$required,parent=if(is.null(s))NA_character_ else a$error$parent,child=if(is.null(s))NA_character_ else a$error$child,required_cells=if(is.null(s))NA_real_ else s$required_cells,available_cells=if(is.null(s))NA_real_ else s$available_cells,supply_ratio=if(is.null(s))NA_real_ else s$supply_ratio,reason=if(is.null(s))NA_character_ else s$passage_failure_reason)
 }))
}
f_ledger<-function(z,candidate,contract){
 do.call(rbind,lapply(names(z$lineages),function(cohort){a<-z$lineages[[cohort]];e<-contract$edges[contract$edges$cohort==cohort&contract$edges$kind=='explicit_parent',,drop=FALSE]
  do.call(rbind,lapply(seq_len(nrow(e)),function(i){r<-if(i<=length(a$edges))a$edges[[i]]else NULL;s<-r$selection;t<-r$target
   data.frame(candidate=candidate,cohort=cohort,parent=e$parent[i],child=e$child[i],status=if(is.null(r))'NOT_REACHED'else if(isTRUE(s$passage_executed))'PASSED'else'FAILED',target_day=if(is.null(r))NA_real_ else t$selected_day,parent_day=if(is.null(r))NA_real_ else s$selected_day,target_cells=if(is.null(r))NA_real_ else t$available_cells,required_cells=if(is.null(r))NA_real_ else s$required_cells,available_cells=if(is.null(r))NA_real_ else s$available_cells,supply_ratio=if(is.null(r))NA_real_ else s$supply_ratio,boundary_scale=if(is.null(r))NA_real_ else s$boundary_scale)
  }))
 }))
}
f_evaluator<-function(c,out,stage,origin,previous=NULL,backend=NULL){
 st<-new.env(parent=emptyenv());st$n<-0L;st$counts<-if(is.null(previous))0L else previous$count;st$last<-NULL;st$times<-list();st$initial<-list();st$archive<-if(is.null(previous))list()else previous$archive;st$anchors<-if(is.null(previous))integer()else previous$anchors
 maxcalls<-switch(stage,tests=128L,diagnostic=1064L,constructor=12000L,validation=2752L);stopifnot(length(maxcalls)==1)
 if(!is.null(previous))stopifnot(identical(previous$contract_sha,c$hash),identical(previous$origin,origin),st$counts>=0,st$counts<=maxcalls)
 log<-file(file.path(out,'evaluations.tsv'),'a');writeLines(paste(c('call','total','event','F','feasible','seconds','key',paste0('u',1:19)),collapse='\t'),log);flush(log)
 failures<-file(file.path(out,'failure_journal.tsv'),'a');failure_header<-FALSE
 emit<-function(x,event,F=NA,ok=NA,seconds=NA) {writeLines(paste(c(st$n,st$counts,event,F,ok,seconds,sf_hash(as.numeric(x)),sprintf('%.17g',x)),collapse='\t'),log);flush(log)}
 limit<-function()stop(structure(list(message='Registered feasibility callback cap exhausted',call=NULL),class=c('feasibility_budget_stop','error','condition')))
 fn<-function(u){
  if(st$counts>=maxcalls)limit();st$n<-st$n+1L;st$counts<-st$counts+1L;emit(u,'START');start<-proc.time()[3]
  if(stage=='constructor'&&st$n<=200L)st$initial[[st$n]]<-as.numeric(u)
  p<-f_decode(c$env,c$d,u);z<-if(is.null(backend))f_physics(c,p)else backend(u);st$last<-z
  if(stage=='constructor'&&z$feasible){
   physical<-which(f_free(c$d)$parameter!='sigma_kary');distinct<-!length(st$archive)||all(vapply(st$archive,function(a)max(abs(a$unit[physical]-u[physical]))>=1e-6,logical(1)))
   if(distinct){
    isanchor<-length(st$anchors)<8L && (!length(st$anchors)||all(vapply(st$archive[st$anchors],function(a)max(abs(a$unit[physical]-u[physical]))>=.02,logical(1))))
    # Archive compact provenance for every distinct feasible point; preserve full
    # states only for first200 physical points or new anchor candidates.
    keep<-length(st$archive)<200L||isanchor;pred<-z$pred;if(!keep)pred$states<-NULL
    record<-list(unit=as.numeric(u),params=p[cg_roles],key=sf_hash(as.numeric(u)),pred=pred,call=st$counts,origin=origin,stage=basename(out),physical=physical,decision=f_short(z),contract_sha=c$hash)
    st$archive[[length(st$archive)+1L]]<-record;if(isanchor)st$anchors<-c(st$anchors,length(st$archive))
    if(keep)sf_rds(record,file.path(out,'members',paste0(record$key,'.rds')))
   }
  }
  if(!is.null(z$lineages)){q<-f_short(z);q$call<-st$counts;write.table(q,failures,sep='\t',row.names=FALSE,col.names=!failure_header,quote=TRUE);failure_header<<-TRUE;flush(failures)}
  seconds<-proc.time()[3]-start;st$times[[st$n]]<-data.frame(call=st$counts,seconds=seconds,feasible=z$feasible);emit(u,'END',sprintf('%.17g',z$F),z$feasible,sprintf('%.9g',seconds))
  if(stage=='constructor'&&length(st$archive)>=200L&&length(st$anchors)>=8L)stop(structure(list(message='Registered feasible bank target reached',call=NULL),class=c('feasibility_bank_complete','error','condition')))
  z$F
 }
 snapshot<-function(extra=list()){
  sf_emit(if(length(st$times))do.call(rbind,st$times)else data.frame(call=integer(),seconds=double(),feasible=logical()),file.path(out,'callback_timings.tsv'))
  c(list(origin=origin,count=st$counts,archive=st$archive,anchors=st$anchors,rng=if(exists('.Random.seed',.GlobalEnv))get('.Random.seed',.GlobalEnv)else NULL,calls=st$n,contract_sha=c$hash),extra)
 }
 list(fn=fn,state=st,snapshot=snapshot,close=function(){close(log);close(failures)},limit=limit)
}
f_save<-function(out,x){sf_rds(x,file.path(out,'state.rds'));keys<-names(x)[vapply(x,function(v)is.atomic(v)&&length(v)==1,logical(1))];sf_emit(data.frame(key=keys,value=vapply(x[keys],as.character,character(1))),file.path(out,'result.tsv'))}
