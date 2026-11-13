# Deterministic, serializable multiscale polling. No sweep/cycle/window cap.
source_directions <- function(active) {
 blocks<-list(resource=c('log10_o2_S0','log10_kappa_O','log10_eta_o2'),growth=c('log10_lam_max','log10_alpha_o2','gamma_growth'),death=c('log10_mu_hp','gamma_mu','log10_O2_crit','n_O'),chromosome=c('log10_p_mis_base','log10_p_misseg','log10_k_o_mis','log10_p_wgd','buffer_smax','log10_buffer_beta','buffer_n_exp'),clearance='log10_k_clear')
 # Preserve the original parameter spelling in the frozen direction set.
 blocks$chromosome[blocks$chromosome=='buffer_n_exp']<-'log10_buffer_n_exp'
 dirs<-list()
 for(n in active)for(sign in c(-1,1)){v<-setNames(rep(0,length(active)),active);v[n]<-sign;dirs[[paste0(n,'_',sign)]]<-v}
 for(block in names(blocks)) {
  a<-intersect(blocks[[block]],active);if(length(a)<2)next
  for(kind in c('same','alternating'))for(sign in c(-1,1)) {
   v<-setNames(rep(0,length(active)),active);v[a]<-sign*if(kind=='same')rep(1,length(a)) else rep(c(1,-1),length.out=length(a));dirs[[paste(block,kind,sign,sep='_')]]<-v
  }
 };dirs
}
machine_new<-function(p,value,lower,upper,fixed=NULL) {
 list(version=1L,best_p=p,best_value=value,lower=lower,upper=upper,fixed=fixed,
  steps=c(.02,.005,.001,.0001),epsilon=1e-10,phase='observe',after_observe='poll_start',
  scale=1L,cycle=1L,cycle_start=value,certifying=FALSE,total_calls=0L,polls=0L,
  observations=0L,visits=setNames(integer(4),as.character(c(.02,.005,.001,.0001))),
  pending=NULL,last_cycle=list(),previous_gain=NA_real_,termination='running')
}
machine_run<-function(st,evaluate,propose,max_callbacks,deadline=Inf,checkpoint=function(s)NULL,poll_done=function(p)NULL) {
 used<-0L;active<-setdiff(names(st$best_p),c(names(st$fixed),'log10_rho_2N','log10_sigma_burden'))
 dirs<-source_directions(active);span<-st$upper-st$lower
 take<-function(p,kind,step=NA_real_,direction='') {
  if(used>=max_callbacks)stop('CONT_CALL_LIMIT')
  if(as.numeric(Sys.time())>=deadline)stop('CONT_TIME_LIMIT')
  value<-evaluate(p,kind,step,direction,st$total_calls+1L)
  stopifnot(length(value)==1L,!is.na(value))
  used<<-used+1L;st$total_calls<<-st$total_calls+1L
  if(value<st$best_value){st$best_value<<-value;st$best_p<<-p}
  value
 }
 reason<-'running'
 tryCatch({
  repeat {
   if(st$phase=='observe') {
    proposal<-propose(st$best_p,st$best_value)
    if(is.null(proposal) || max(abs(proposal$parameters-st$best_p))<1e-14)st$phase<-st$after_observe
    else {st$pending_observation<-proposal;st$phase<-'observation_evaluate';checkpoint(st)}
   } else if(st$phase=='observation_evaluate') {
    proposal<-st$pending_observation;value<-take(proposal$parameters,'observation')
    stopifnot(abs(value-proposal$expected)<1e-8)
    st$observations<-st$observations+1L;st$pending_observation<-NULL;st$phase<-st$after_observe;checkpoint(st)
   } else if(st$phase=='poll_start') {
    base<-st$best_p;step<-st$steps[st$scale];candidates<-list()
    for(name in names(dirs)) {
     q<-base;q[active]<-base[active]+step*span[active]*dirs[[name]]
     if(any(q<st$lower | q>st$upper))next
     candidates[[length(candidates)+1L]]<-list(parameters=q,direction=name)
    }
    st$pending<-list(center=base,before=st$best_value,step=step,candidates=candidates,next_index=1L,evaluated=list(),cycle=st$cycle,certifying=st$certifying)
    st$phase<-'poll_evaluate'
   } else if(st$phase=='poll_evaluate') {
    while(st$pending$next_index<=length(st$pending$candidates)) {
     j<-st$pending$next_index;q<-st$pending$candidates[[j]]
     value<-take(q$parameters,'poll',st$pending$step,q$direction)
     st$pending$evaluated[[j]]<-c(q,list(value=value,sequence=st$total_calls))
     st$pending$next_index<-j+1L
    }
    st$polls<-st$polls+1L;st$visits[st$scale]<-st$visits[st$scale]+1L
    poll_done(c(st$pending,list(poll=st$polls,after=st$best_value)))
    st$last_cycle[[st$scale]]<-list(center=st$pending$center,before=st$pending$before,after=st$best_value,step=st$pending$step,poll=st$polls)
    improved<-st$best_value<st$pending$before
    st$pending<-NULL;st$after_observe<-'advance';st$phase<-if(improved)'observe' else 'advance';checkpoint(st)
   } else if(st$phase=='advance') {
    if(st$scale<length(st$steps)){st$scale<-st$scale+1L;st$phase<-'poll_start'}
    else {
     gain<-st$cycle_start-st$best_value;st$previous_gain<-gain
     if(st$certifying && gain<=st$epsilon){reason<-'poll_stationary';st$termination<-reason;checkpoint(st);break}
     st$certifying<-gain<=st$epsilon;st$cycle<-st$cycle+1L;st$scale<-1L;st$cycle_start<-st$best_value;st$last_cycle<-list();st$phase<-'poll_start';checkpoint(st)
    }
   } else stop('Unknown checkpoint phase')
  }
 },error=function(e) {
  if(grepl('CONT_CALL_LIMIT',conditionMessage(e),fixed=TRUE))reason<<-'call_limit'
  else if(grepl('CONT_TIME_LIMIT',conditionMessage(e),fixed=TRUE)||grepl('IVI_SOFT_TIME_LIMIT',conditionMessage(e),fixed=TRUE))reason<<-'time_limit'
  else stop(e)
 })
 st$termination<-reason;checkpoint(st)
 list(state=st,reason=reason,window_calls=used,improved=st$cycle_start-st$best_value>st$epsilon || (!is.na(st$previous_gain)&&st$previous_gain>st$epsilon))
}
