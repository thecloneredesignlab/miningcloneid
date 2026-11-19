# One source-bound finite allocation; no observation loss or fitting dispatch.
main<-function(){
 args<-commandArgs(TRUE);study<-normalizePath(args[1]);task_id<-args[2];task<-sub('_R1$','',task_id);out<-file.path(study,'verification',task_id)
 for(name in c('production_bindings.R','interface.R','core.R','setup.R','tests.R','diagnose.R','build.R','validate.R'))source(file.path(study,'release',name))
 on.exit(sf_inventory(out),add=TRUE)
 checks<-list();check<-function(n,ok,d=''){checks[[length(checks)+1L]]<<-data.frame(check=n,passed=isTRUE(ok),detail=d);sf_emit(do.call(rbind,checks),file.path(out,'checks.tsv'));if(!isTRUE(ok))stop('Check failed: ',n)}
 env<-sf_backend(study,out,native=!task%in%c('P0','P1C'))
 if(task=='P0')f_setup(study,out,env,check)else if(task=='P1')f_tests(study,out,env,check)else if(task=='P1C')f_canary(study,out,env,check)else if(startsWith(task,'P2_'))f_diagnose(study,task,out,env,check)else if(grepl('_(check1|check2|robust)$',task))f_validate(study,task,out,env,check)else f_build(study,task,out,env,check)
 sf_emit(data.frame(status='PASS',task=task_id,D19_fit=0,N20=0,profiles=0,tumor=0),file.path(out,'application_result.tsv'));sf_emit(data.frame(status='PASS',task=task_id),file.path(out,'gate.tsv'))
}
tryCatch(main(),error=function(e){cat('FEASIBILITY_PHASE_FAILED: ',conditionMessage(e),'\n',sep='');quit(status=1,save='no')})
