# Registered blind population-linked recovery stages and exact-vector checks.
r_fit<-function(study,task,out,env,check){
 origin<-sub('_[^_]+$','',task);stage<-sub('^.*_','',task);kind<-substr(origin,1,1)
 case<-readRDS(file.path(study,'verification/I01/cases',if(kind=='D')'D19.rds' else 'N20.rds'));ctx<-r_context(env,case)
 chain<-c('screen','bench','de1b',paste0('de',2:6),'local','pattern','verify1','verify2');idx<-match(stage,chain);stopifnot(!is.na(idx))
 previous<-if(idx>1)readRDS(file.path(study,'verification',paste0(origin,'_',chain[idx-1]),'state.rds'))else NULL
 bank<-readRDS(file.path(study,'verification/I01/banks',paste0(origin,'.rds')))
 assign('.Random.seed',if(is.null(previous))bank$state else previous$rng,.GlobalEnv)
 ev<-r_evaluator(ctx,out,origin,stage,previous);on.exit(ev$close(),add=TRUE)
 completed<-FALSE;on.exit(if(!completed)r_result(out,ev$snapshot(list(status='FAILED_PARTIAL_ATTEMPTS_CHARGED'))),add=TRUE)
 finish<-function(extra){r_result(out,ev$snapshot(extra));completed<<-TRUE}
 same_pop<-function(a,b)identical(sort(apply(a,1,sf_hash)),sort(apply(b,1,sf_hash)))
 if(stage=='screen'){
  raw<-list();proposals<-list();feasible<-integer()
  for(i in 1:1000){z<-r_proposal(ctx);raw[[i]]<-z$unit;proposals[[i]]<-data.frame(draw=i,proposal=z$kind)
   ev$fn(z$unit);if(ev$state$last$feasible)feasible<-c(feasible,i);if(length(feasible)>=8)break}
  if(length(feasible)<8){finish(list(status='SCIENTIFIC_SCREEN_INSUFFICIENT',feasible=length(feasible)));return(invisible(NULL))}
  while(length(raw)<200){z<-r_proposal(ctx);raw[[length(raw)+1L]]<-z$unit;proposals[[length(raw)]]<-data.frame(draw=length(raw),proposal=paste0(z$kind,'_unscored_fill'))}
  pop<-do.call(rbind,raw[1:200]);pop[1:8,]<-do.call(rbind,raw[feasible[1:8]])
  sf_rds(list(raw=raw,first8=feasible[1:8],initialpop=pop),file.path(out,'raw_bank.rds'));sf_emit(do.call(rbind,proposals),file.path(out,'proposals.tsv'))
  sf_emit(data.frame(row=1:8,source_draw=feasible[1:8],duplicate_existing=vapply(feasible[1:8],function(i)i<=200,logical(1))),file.path(out,'population_replacements.tsv'))
  finish(list(status='READY',population=pop,feasible=length(feasible)));return(invisible(NULL))
 }
 if(stage=='bench'||startsWith(stage,'de')){
  suppressPackageStartupMessages(library(DEoptim));pop<-previous$population;generations<-if(stage=='bench')1L else if(stage=='de1b')59L else 60L
  stopifnot(is.matrix(pop),identical(dim(pop),c(200L,nrow(ctx$free))))
  sf_rds(list(population=pop,rng=.Random.seed),file.path(out,'requested_population.rds'))
  de<-DEoptim::DEoptim(ev$fn,rep(0,nrow(ctx$free)),rep(1,nrow(ctx$free)),control=list(NP=200L,itermax=generations,trace=FALSE,initialpop=pop,strategy=2L,F=.8,CR=.9,reltol=0,steptol=generations+1L,VTR=-Inf,parallelType=0))
  actual<-do.call(rbind,ev$state$initial);sf_rds(actual,file.path(out,'actual_initial_population.rds'))
  check('actual_requested_native_population',same_pop(pop,actual));check('native_callbacks_reconcile',de$optim$nfeval==ev$state$n&&ev$state$n==200L*(generations+1L))
  sf_rds(de,file.path(out,'native_DE_result.rds'));finish(list(status='READY',population=de$member$pop,native_evaluations=de$optim$nfeval,generations=generations));return(invisible(NULL))
 }
 x<-previous$best$unit;stopifnot(is.numeric(x),length(x)==nrow(ctx$free),is.finite(previous$best$score))
 replay<-ev$fn(x);check('native_vector_score_handoff',abs(replay-previous$best$score)<=1e-8)
 check('native_vector_parameter_handoff',identical(r_decode(ctx,x)[cg_roles],previous$best$params))
 if(stage=='local'){
  suppressPackageStartupMessages(library(nloptr));stopped<-FALSE
  nl<-tryCatch(nloptr::nloptr(x,ev$fn,lb=rep(0,length(x)),ub=rep(1,length(x)),opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=6000L,maxtime=0,xtol_abs=rep(1e-5,length(x)),ftol_abs=if(kind=='D')1e-8 else 1e-6,print_level=0)),recovery_budget_stop=function(e){stopped<<-TRUE;NULL})
  sf_rds(nl,file.path(out,'native_Subplex_result.rds'));finish(list(status='READY',wrapper_allowance_exhausted=stopped,native_status=if(is.null(nl))NA_integer_ else nl$status));return(invisible(NULL))
 }
 if(stage=='pattern'){
  pol<-r_poll(ev$fn,x,r_directions(ctx),c(.02,.005,.001,.0001),1e-8,min(1999L,2000L-ev$state$counts['pattern']))
  sf_emit(pol$rows,file.path(out,'pattern_polls.tsv'));finish(list(status='READY',pattern_budget_limited=pol$budget_limited));return(invisible(NULL))
 }
 # Verification is at one frozen exact vector. Improvements are saved but do not
 # inherit this vector's replay/poll certificate. A second pass can check them.
 anchor<-list(score=replay,unit=x,key=sf_hash(x),params=r_decode(ctx,x)[cg_roles],pred=ev$state$last$pred,loss=ev$state$last$loss)
 original<-ev$state$last;again<-ev$fn(x);fresh<-ev$state$last
 check('exact_vector_repeat_score',abs(replay-again)<=1e-8)
 check('exact_protocol_decisions',identical(original$pred$selection,fresh$pred$selection))
 check('exact_vector_state_repeat',all(vapply(names(original$pred$states),function(n)all(abs(original$pred$states[[n]]-fresh$pred$states[[n]])<=1e-10*pmax(1,abs(original$pred$states[[n]]),abs(fresh$pred$states[[n]]))),logical(1))))
 check('repeat_probability_TV',max(vapply(names(original$pred$kary),function(n)sum(abs(original$pred$kary[[n]]-fresh$pred$kary[[n]]))/2,numeric(1)))<=1e-10)
 if(kind=='N'){full<-cg_score(env,case$design,fresh$runs,r_decode(ctx,x),fresh$pred,case$data);check('independent_original_score',abs(full$objective-replay)<=1e-8)}else{
  g<-fresh$pred$growth[names(case$data$growth)]-case$data$growth
  d<-sum(g*g)/(2*length(g)*case$known_sigma_growth^2)+sum(vapply(names(case$data$prob),function(n){a<-fresh$pred$kary[[n]];a<-a/sum(a);sum((sqrt(a)-sqrt(case$data$prob[[n]]))^2)},numeric(1)))/(2*length(case$data$prob))
  check('independent_expected_distance',abs(d-replay)<=1e-8)
 }
 directions<-r_directions(ctx);rows<-list();seen<-character();complete<-TRUE;max_gain<-0
 for(h in c(.001,.0001,.00002))for(j in seq_len(nrow(directions)))for(sign in c(-1,1)){
  u<-r_admissible(x,sign*directions[j,],h);if(is.null(u))next;k<-sf_hash(u);if(k%in%seen)next;seen<-c(seen,k)
  q<-tryCatch(ev$fn(u),recovery_budget_stop=function(e)NA_real_)
  if(is.na(q)){complete<-FALSE;next};max_gain<-max(max_gain,replay-q)
  rows[[length(rows)+1L]]<-data.frame(mesh=h,direction=j,sign=sign,key=k,score=q,improvement=replay-q)
 }
 sigma_gain<-0
 if(kind=='N'){
  residual<-case$data$growth-anchor$pred$growth[names(case$data$growth)];tab<-ctx$free;ii<-match('sigma_growth',tab$parameter);sc<-min(tab$upper[ii],max(tab$lower[ii],sqrt(mean(residual^2))))
  u<-x;u[ii]<-(env$o2prof_transform(sc,tab[ii,,drop=FALSE])-ctx$low[ii])/(ctx$high[ii]-ctx$low[ii]);u[ii]<-max(0,min(1,u[ii]))
  analytic<-r_loss(ctx,anchor$pred,modifyList(r_decode(ctx,x),list(sigma_growth=sc)))$value;sigma_gain<-replay-analytic
  actual<-tryCatch(ev$fn(u),recovery_budget_stop=function(e)NA_real_)
  if(is.na(actual))complete<-FALSE else check('single_vector_variance_check',abs(actual-analytic)<=1e-8)
  sf_emit(data.frame(anchor_key=anchor$key,sigma_cond=sc,analytic_score=analytic,full_score=actual,improvement=sigma_gain),file.path(out,'variance_check.tsv'))
 }
 sf_emit(if(length(rows))do.call(rbind,rows)else data.frame(mesh=numeric(),direction=integer(),sign=integer(),key=character(),score=numeric(),improvement=numeric()),file.path(out,'verification_polls.tsv'))
 predpass<-kind=='D'&&anchor$score<=1e-4&&anchor$loss$rms_growth<=.01&&anchor$loss$max_growth<=.03&&anchor$loss$max_kary_tv<=.005
 localpass<-complete&&max_gain<=1e-5&&sigma_gain<=1e-5
 sf_rds(anchor,file.path(out,'verified_anchor.rds'))
 sf_emit(data.frame(origin=origin,case=kind,key=anchor$key,score=anchor$score,best_found=ev$state$best$score,best_found_key=ev$state$best$key,verified=TRUE,complete_poll=complete,local_pass=localpass,maximum_gain=max_gain,sigma_gain=sigma_gain,prediction_pass=predpass,rms_growth=anchor$loss$rms_growth,max_growth=anchor$loss$max_growth,max_kary_tv=anchor$loss$max_kary_tv),file.path(out,'quality.tsv'))
 finish(list(status='VERIFIED',verified_anchor=anchor,local_pass=localpass,prediction_pass=predpass,complete_poll=complete))
}
