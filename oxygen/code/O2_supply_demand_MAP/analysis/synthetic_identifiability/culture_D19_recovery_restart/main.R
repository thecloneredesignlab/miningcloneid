main<-function(){
 args<-commandArgs(TRUE);study<-normalizePath(args[1]);task<-args[2];out<-file.path(study,'verification',task)
 for(name in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R','core.R','engineering.R','fit.R','build.R','validate.R'))source(file.path(study,'release',name))
 on.exit(sf_inventory(out),add=TRUE)
 checks<-list();check<-function(n,ok,d=''){checks[[length(checks)+1L]]<<-data.frame(check=n,passed=isTRUE(ok),detail=d);sf_emit(do.call(rbind,checks),file.path(out,'checks.tsv'));if(!isTRUE(ok))stop('Check failed: ',n)}
 env<-sf_backend(study,out,native=!task%in%c('T00','T02D','T03F'))
 if(task=='T00')d_setup(study,out,env,check)else if(task=='T01')d_tests(study,out,env,check)else if(task%in%c('T02D','T03F'))d_canaries(study,out,env,check,if(task=='T02D')'D'else'F')else if(startsWith(task,'D_'))d_fit(study,task,out,env,check)else if(grepl('_(check1|check2|robust)$',task))f_validate(study,task,out,env,check)else if(startsWith(task,'F_'))f_build(study,task,out,env,check)else stop('Unauthorized task')
 sf_emit(data.frame(status='PASS',task=task,N20=0,profiles=0,tumor=0),file.path(out,'application_result.tsv'));sf_emit(data.frame(status='PASS',task=task),file.path(out,'gate.tsv'))
}
tryCatch(main(),error=function(e){cat('D19_PHASE_FAILED: ',conditionMessage(e),'\n',sep='');quit(status=1,save='no')})
