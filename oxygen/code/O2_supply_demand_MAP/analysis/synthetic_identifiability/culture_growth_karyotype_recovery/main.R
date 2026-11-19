# One isolated, manifest-bound phase. Unrecognized errors are never penalties.
main<-function(){
 args<-commandArgs(TRUE);study<-normalizePath(args[1]);task_id<-args[2];task<-sub('_R1$','',task_id);out<-file.path(study,'verification',task_id)
 for(f in c('production_bindings.R','interface.R','core.R','engineering.R','fit.R'))source(file.path(study,'release',f))
 on.exit(sf_inventory(out),add=TRUE)
 env<-sf_backend(study,out,native=task=='I02'||grepl('^[DN]_',task))
 checks<-list();check<-function(n,ok,d=''){checks[[length(checks)+1L]]<<-data.frame(check=n,passed=isTRUE(ok),detail=d);sf_emit(do.call(rbind,checks),file.path(out,'checks.tsv'));if(!isTRUE(ok))stop('Check failed: ',n)}
 if(task%in%c('I01','I02','I03'))r_engineering(study,task,out,env,check)else r_fit(study,task,out,env,check)
 sf_emit(data.frame(status='PASS',task=task),file.path(out,'application_result.tsv'));sf_emit(data.frame(status='PASS',task=task),file.path(out,'gate.tsv'))
}
tryCatch(main(),error=function(e){cat('RECOVERY_PHASE_FAILED: ',conditionMessage(e),'\n',sep='');quit(status=1,save='no')})
