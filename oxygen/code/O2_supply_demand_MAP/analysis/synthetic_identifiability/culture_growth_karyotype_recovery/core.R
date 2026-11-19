# Finite candidate evaluator. No private truth access is needed by fitting workers.
r_limits<-c(screen=1000L,global=73200L,local=6000L,pattern=2000L,verify=800L,handoff=1000L,external=512L)
r_context<-function(env,case){
 t<-case$design$bounds;free<-t[if(case$kind=='D')t$parameter!='sigma_growth' else rep(TRUE,nrow(t)),,drop=FALSE]
 low<-env$o2prof_transform(free$lower,free);high<-env$o2prof_transform(free$upper,free)
 h<-if(case$kind=='N')lapply(case$data$kary,function(v)tabulate(match(v,22:154),133))else NULL
 list(env=env,case=case,free=free,low=low,high=high,core=env$build_model_core(cfg=case$design$cfg),hist=h)
}
r_decode<-function(c,u){
 stopifnot(is.numeric(u),length(u)==nrow(c$free),all(is.finite(u)),all(u>=0),all(u<=1))
 t<-c$case$design$bounds;p<-cg_neutral(c$env,t);if(c$case$kind=='D')p$sigma_growth<-c$case$known_sigma_growth
 z<-setNames(vapply(seq_len(nrow(t)),function(i)c$env$o2prof_transform(p[[t$parameter[i]]],t[i,,drop=FALSE]),numeric(1)),t$coordinate)
 z[c$free$coordinate]<-c$low+as.numeric(u)*(c$high-c$low)
 cg_from_free(c$env,c$case$design,z)
}
r_encode<-function(c,p){
 t<-c$free;z<-vapply(seq_len(nrow(t)),function(i)c$env$o2prof_transform(p[[t$parameter[i]]],t[i,,drop=FALSE]),numeric(1))
 as.numeric((z-c$low)/(c$high-c$low))
}
r_loss<-function(c,pred,p){
 if(c$case$kind=='D'){
  g<-(pred$growth[names(c$case$data$growth)]-c$case$data$growth)/c$case$known_sigma_growth
  tv<-vapply(names(c$case$data$prob),function(id)sum(abs(pred$kary[[id]]/sum(pred$kary[[id]])-c$case$data$prob[[id]]))/2,numeric(1))
  hk<-vapply(names(c$case$data$prob),function(id)sum((sqrt(pred$kary[[id]]/sum(pred$kary[[id]]))-sqrt(c$case$data$prob[[id]]))^2)/2,numeric(1))
  list(value=mean(g*g)/2+mean(hk),growth=mean(g*g)/2,kary=mean(hk),rms_growth=sqrt(mean(g*g)),max_growth=max(abs(g)),max_kary_tv=max(tv))
 }else{
  residual<-c$case$data$growth-pred$growth[names(c$case$data$growth)]
  g<-mean(-log(p$sigma_growth)-log(2*pi)/2-residual^2/(2*p$sigma_growth^2))
  k<-mean(vapply(names(c$hist),function(id)sum(c$hist[[id]]*log(pmax(pred$kary[[id]],1e-12)))/sum(c$hist[[id]]),numeric(1)))
  list(value=-(g+k),growth=g,kary=k,rms_growth=NA_real_,max_growth=NA_real_,max_kary_tv=NA_real_)
 }
}
r_physics<-function(c,p){
 runs<-tryCatch(lapply(c$case$design$adapters,function(a)c$env$ivt_run_lineage(a,cfg=c$case$design$cfg,run_params=p,model_core=c$core)),invitro_protocol_infeasible=function(e)NULL)
 if(is.null(runs))return(list(feasible=FALSE))
 pred<-cg_predict(c$env,c$case$design,runs,p);loss<-r_loss(c,pred,p)
 stopifnot(is.finite(loss$value));list(feasible=TRUE,runs=runs,pred=pred,loss=loss)
}
r_activity<-function(p){
 inactive<-character()
 if(p$p_misseg==0)inactive<-c(inactive,'k_o_mis')
 if(p$buffer_beta==0)inactive<-c(inactive,'buffer_n_exp')
 if(p$buffer_smax==0)inactive<-c(inactive,'buffer_beta','buffer_n_exp')
 unique(inactive)
}
r_directions<-function(c){
 n<-nrow(c$free);D<-diag(n);rownames(D)<-c$free$parameter
 pairs<-list(c('lam_max','alpha_o2',1),c('mu_hp','k_o_mis',1),c('p_misseg','p_mis_base',-1),c('buffer_smax','buffer_beta',1),c('init_sd_2N','sigma_kary',-1),c('init_sd_4N','sigma_kary',-1))
 for(pair in pairs){v<-numeric(n);v[match(pair[1],c$free$parameter)]<-1/sqrt(2);v[match(pair[2],c$free$parameter)]<-as.numeric(pair[3])/sqrt(2);D<-rbind(D,v)}
 colnames(D)<-c$free$parameter;D
}
r_admissible<-function(x,d,h){
 dist<-c((1-x[d>0])/d[d>0],-x[d<0]/d[d<0]);step<-min(c(h,dist))
 if(!is.finite(step)||step<=0)return(NULL)
 y<-x+step*d;tol<-64*.Machine$double.eps
 if(any(y < -tol|y>1+tol))stop('Material boundary error')
 y<-pmin(1,pmax(0,y));if(identical(as.numeric(y),as.numeric(x)))NULL else as.numeric(y)
}
r_evaluator<-function(c,out,origin,stage,previous=NULL){
 state<-new.env(parent=emptyenv());state$counts<-setNames(integer(length(r_limits)),names(r_limits));state$n<-0L;state$last<-NULL;state$best<-NULL;state$timings<-list();state$initial<-list()
 if(!is.null(previous)){state$counts<-previous$counts;state$best<-previous$best}
 stopifnot(sum(state$counts[setdiff(names(r_limits),'external')])<=84000)
 if(is.null(state$best))state$best<-list(score=Inf)
 journal<-file(file.path(out,'evaluations.tsv'),'a');header<-c('call','event','bucket','score','seconds','feasible','candidate_sha',paste0('u',seq_len(nrow(c$free))))
 writeLines(paste(header,collapse='\t'),journal);flush(journal)
 primary<-switch(stage,screen='screen',bench='global',de1b='global',de2='global',de3='global',de4='global',de5='global',de6='global',local='local',pattern='pattern',verify1='verify',verify2='verify',external='external')
 stopifnot(length(primary)==1)
 capstop<-function()stop(structure(list(message='Registered candidate allowance exhausted',call=NULL),class=c('recovery_budget_stop','error','condition')))
 fn<-function(u){
  bucket<-if(stage=='de1b'&&state$n<200L)'handoff'else primary
  if(state$counts[bucket]>=r_limits[bucket] || (stage!='external'&&sum(state$counts[setdiff(names(r_limits),'external')])>=84000))capstop()
  state$n<-state$n+1L;state$counts[bucket]<-state$counts[bucket]+1L;start<-proc.time()[3]
  p<-r_decode(c,u);key<-sf_hash(as.numeric(u));line<-c(state$n,'START',bucket,NA,NA,NA,key,sprintf('%.17g',u));writeLines(paste(line,collapse='\t'),journal);flush(journal)
  if((stage=='bench'||startsWith(stage,'de'))&&state$n<=200L)state$initial[[length(state$initial)+1L]]<-as.numeric(u)
  z<-r_physics(c,p);state$last<-z;score<-if(z$feasible)z$loss$value else Inf
  if(is.finite(score)&&score<state$best$score){
   pred<-z$pred;if(!startsWith(stage,'verify'))pred$states<-NULL;pred$latent_kary<-NULL
   state$best<-list(score=score,unit=as.numeric(u),params=p[cg_roles],pred=pred,loss=z$loss,key=key,origin=origin,source_stage=stage,data_sha=c$case$data_sha)
   sf_rds(state$best,file.path(out,'best.rds'));dir.create(file.path(out,'improvements'),showWarnings=FALSE)
   sf_rds(state$best,file.path(out,'improvements',sprintf('%06d.rds',state$n)))
  }
  seconds<-proc.time()[3]-start;state$timings[[state$n]]<-data.frame(call=state$n,seconds=seconds,feasible=z$feasible)
  line<-c(state$n,'END',bucket,sprintf('%.17g',score),sprintf('%.9g',seconds),z$feasible,key,sprintf('%.17g',u));writeLines(paste(line,collapse='\t'),journal);flush(journal)
  score
 }
 snapshot<-function(extra=list()){
  sf_emit(if(length(state$timings))do.call(rbind,state$timings)else data.frame(call=integer(),seconds=double(),feasible=logical()),file.path(out,'callback_timings.tsv'))
  c(list(origin=origin,stage=stage,counts=state$counts,best=state$best,rng=if(exists('.Random.seed',.GlobalEnv))get('.Random.seed',.GlobalEnv)else NULL,calls=state$n,data_sha=c$case$data_sha),extra)
 }
 list(fn=fn,state=state,snapshot=snapshot,close=function()close(journal),capstop=capstop)
}
r_poll<-function(fn,x,directions,meshes,tolerance,max_calls){
 calls<-0L;q<-fn(x);calls<-calls+1L;rows<-list();exhausted<-FALSE
 for(h in meshes){repeat{
  oldq<-q;winner<-x;winq<-q;seen<-character();complete<-TRUE;maximum<-0
  for(j in seq_len(nrow(directions)))for(sign in c(-1,1)){
   y<-r_admissible(x,sign*directions[j,],h);if(is.null(y))next;key<-sf_hash(y);if(key%in%seen)next;seen<-c(seen,key)
   if(calls>=max_calls){complete<-FALSE;exhausted<-TRUE;break}
   z<-tryCatch(fn(y),recovery_budget_stop=function(e){exhausted<<-TRUE;NA_real_});if(is.na(z)){complete<-FALSE;break};calls<-calls+1L
   maximum<-max(maximum,q-z);if(z<winq){winner<-y;winq<-z}
  }
  x<-winner;q<-winq;rows[[length(rows)+1L]]<-data.frame(mesh=h,complete=complete,before=oldq,after=q,maximum_improvement=maximum,calls=calls)
  if(exhausted||oldq-q<=tolerance)break
 }
 if(exhausted)break}
 list(unit=x,score=q,rows=do.call(rbind,rows),budget_limited=exhausted,calls=calls)
}
r_proposal<-function(c){
 x<-runif(nrow(c$free));kind<-'uniform'
 if(runif(1)<.5){kind<-'feasibility_mixture';x[match('lam_max',c$free$parameter)]<-runif(1,.5,1);x[match('mu_hp',c$free$parameter)]<-runif(1,0,.5)}
 list(unit=x,kind=kind)
}
r_result<-function(out,x){sf_rds(x,file.path(out,'state.rds'));sf_emit(data.frame(key=names(x)[vapply(x,function(v)is.atomic(v)&&length(v)==1,logical(1))],value=vapply(x[vapply(x,function(v)is.atomic(v)&&length(v)==1,logical(1))],as.character,character(1))),file.path(out,'result.tsv'))}
