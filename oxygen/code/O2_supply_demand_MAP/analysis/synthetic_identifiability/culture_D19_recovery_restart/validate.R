# Fresh unmodified production replays and nonselective local-feasibility probes.
f_validate<-function(study,task,out,env,check){
 origin<-sub('_(check1|check2|robust)$','',task);stage<-sub('^.*_','',task)
 path<-readLines(file.path(study,'controller',paste0(origin,'_bank_path.txt')));stopifnot(length(path)==1,startsWith(path,paste0(origin,'_')))
 bank<-readRDS(file.path(study,'verification',path,'bank.rds'));contract<-readRDS(file.path(study,'inputs/public/contract.rds'));c<-f_context(env,contract)
 check('bank_contract_identity',identical(bank$contract_sha,sf_hash(contract))&&bank$origin==origin&&length(bank$members)==200)
 check('physical_unique_members',!anyDuplicated(vapply(bank$members,`[[`,character(1),'key')))
 # Per-family validation follows exactly200+200+at most288; aggregate<=2752.
 previous<-if(stage=='check1')NULL else readRDS(file.path(study,'verification',paste0(origin,'_',if(stage=='check2')'check1'else'check2'),'state.rds'))
 if(!is.null(previous)){previous$count<-previous$validation_count;previous$archive<-list();previous$anchors<-integer()}
 backend<-function(u){p<-f_decode(env,c$d,u);r<-f_reference(c,p);r$F<-if(r$feasible)0 else 1.1;r}
 ev<-f_evaluator(c,out,'validation',origin,previous,backend);on.exit(ev$close(),add=TRUE);complete<-FALSE
 on.exit(if(!complete)f_save(out,ev$snapshot(list(status='FAILED_PARTIAL',validation_count=ev$state$counts))),add=TRUE)
 rows<-list()
 if(stage%in%c('check1','check2')){
  physical<-which(f_free(c$d)$parameter!='sigma_kary');u<-do.call(rbind,lapply(bank$members,`[[`,'unit'))
  check('distinct200_physical',all(vapply(1:199,function(i)all(apply(abs(sweep(u[(i+1):200,physical,drop=FALSE],2,u[i,physical],'-')),1,max)>=1e-6),logical(1))))
  check('eight_physical_anchors',all(vapply(1:7,function(i)all(apply(abs(sweep(u[(i+1):8,physical,drop=FALSE],2,u[i,physical],'-')),1,max)>=.02),logical(1))))
  for(i in 1:200){m<-bank$members[[i]];ev$fn(m$unit);r<-ev$state$last;check(paste0('complete_feasible_',i),r$feasible);f_compare(m$pred,r$pred,check,paste0('member_',i));rows[[i]]<-data.frame(row=i,key=m$key,status='PASS')}
 }else{
  physical<-which(f_free(c$d)$parameter!='sigma_kary')
  for(i in 1:8)for(j in physical)for(sign in c(-1,1)){
   u<-bank$members[[i]]$unit;h<-min(.0001,if(sign>0)1-u[j]else u[j]);if(h<=0)next;u[j]<-u[j]+sign*h;if(u[j]<0&&u[j]>=-64*.Machine$double.eps)u[j]<-0;if(u[j]>1&&u[j]<=1+64*.Machine$double.eps)u[j]<-1
   ev$fn(u);rows[[length(rows)+1L]]<-data.frame(anchor=i,parameter=f_free(c$d)$parameter[j],sign=sign,actual_step=h,feasible=ev$state$last$feasible,key=sf_hash(as.numeric(u)))
  }
 }
 check('per_family_validation_cap',ev$state$counts<=688)
 sf_emit(do.call(rbind,rows),file.path(out,if(stage=='robust')'bank_robustness.tsv'else'bank_replay_checks.tsv'))
 f_save(out,ev$snapshot(list(status='PASS',validation_count=ev$state$counts)));complete<-TRUE
}
