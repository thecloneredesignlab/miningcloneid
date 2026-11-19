n_candidate<-function(c,u,z,origin,stage){
 p<-n_decode(c,u);pred<-z$pred;pred$states<-NULL
 list(score=z$loss$value,unit=as.numeric(u),params=p[cg_roles],pred=pred,loss=z$loss,key=sf_hash(as.numeric(u)),origin=origin,source_stage=stage,identity=c$identity,data_sha=c$case$data_sha)
}
n_validate_state<-function(c,x,origin){
 stopifnot(x$origin==origin,identical(x$identity,c$identity),identical(x$data_sha,c$case$data_sha),identical(x$counts,setNames(as.integer(x$counts),names(n_limits))),all(x$counts>=0),all(x$counts<=n_limits),sum(x$counts)<=84000,identical(sf_hash(x$rng),x$rng_sha))
 if(is.finite(x$best$score))stopifnot(identical(x$best$params,n_decode(c,x$best$unit)[cg_roles]),identical(x$best$key,sf_hash(x$best$unit)),x$best$identity==c$identity)
 invisible(TRUE)
}
n_evaluator<-function(c,out,origin,stage,previous=NULL,backend=NULL){
 st<-new.env(parent=emptyenv());st$counts<-setNames(integer(length(n_limits)),names(n_limits));st$n<-0L;st$best<-list(score=Inf);st$last<-NULL;st$times<-list();st$initial<-list()
 if(!is.null(previous)){n_validate_state(c,previous,origin);st$counts<-previous$counts;st$best<-previous$best}
 primary<-if(stage=='bench'||startsWith(stage,'de'))'global'else if(stage%in%c('local','pattern','joint','restore'))stage else 'verify'
 log<-file(file.path(out,'evaluations.tsv'),'a');writeLines(paste(c('call','event','bucket','score','seconds','feasible','key',paste0('u',1:20)),collapse='\t'),log);flush(log)
 fail<-file(file.path(out,'failure_journal.tsv'),'a');failure_header<-FALSE
 capstop<-function()stop(structure(list(message='Registered N20 callback allowance exhausted',call=NULL),class=c('recovery_budget_stop','error','condition')))
 fn<-function(u,bucket=NULL){
  if(is.null(bucket))bucket<-primary
  stopifnot(bucket%in%names(n_limits));if(st$counts[bucket]>=n_limits[bucket]||sum(st$counts)>=84000)capstop()
  st$n<-st$n+1L;st$counts[bucket]<-st$counts[bucket]+1L;t<-proc.time()[3];key<-sf_hash(as.numeric(u))
  emit<-function(event,score=NA,seconds=NA,ok=NA){writeLines(paste(c(st$n,event,bucket,sprintf('%.17g',score),sprintf('%.9g',seconds),ok,key,sprintf('%.17g',u)),collapse='\t'),log);flush(log)}
  emit('START');p<-n_decode(c,u)
  if((stage=='bench'||startsWith(stage,'de'))&&st$n<=200L)st$initial[[st$n]]<-as.numeric(u)
  z<-if(is.null(backend))n_physics(c,p,reference=endsWith(stage,'b')&&startsWith(stage,'replay'))else backend(u);st$last<-z
  score<-if(z$feasible)z$loss$value else Inf;stopifnot(!is.na(score),score>=0)
  if(!is.null(z$lineages)){
   q<-f_short(z);q$call<-st$n;write.table(q,fail,sep='\t',row.names=FALSE,col.names=!failure_header,quote=TRUE);failure_header<<-TRUE;flush(fail)
   if(!z$feasible){sf_rds(lapply(z$lineages,function(a)a$error),file.path(out,'failures',sprintf('%06d.rds',st$n)));sf_emit(f_ledger(z,key,c$contract),file.path(out,'failures',sprintf('%06d_edges.tsv',st$n)))}
  }
  if(is.finite(score)){
   improves<-score<st$best$score
   if(improves){
    candidate<-n_candidate(c,u,z,origin,stage)
    st$best<-candidate
    sf_rds(candidate,file.path(out,'best.rds'))
    sf_rds(candidate,file.path(out,'improvements',sprintf('%06d.rds',st$n)))
   }
  }
  seconds<-proc.time()[3]-t;st$times[[st$n]]<-data.frame(call=st$n,seconds=seconds,feasible=z$feasible,score=score);emit('END',score,seconds,z$feasible);score
 }
 snapshot<-function(extra=list()){
  sf_emit(if(length(st$times))do.call(rbind,st$times)else data.frame(call=integer(),seconds=double(),feasible=logical(),score=double()),file.path(out,'callback_timings.tsv'))
  rng<-get('.Random.seed',.GlobalEnv)
  c(list(origin=origin,stage=stage,counts=st$counts,best=st$best,rng=rng,rng_sha=sf_hash(rng),calls=st$n,identity=c$identity,data_sha=c$case$data_sha),extra)
 }
 list(fn=fn,state=st,snapshot=snapshot,close=function(){close(log);close(fail)},capstop=capstop)
}
n_bank<-function(study,origin){
 readRDS(file.path(study,'inputs/banks',paste0(origin,'.rds')))
}
