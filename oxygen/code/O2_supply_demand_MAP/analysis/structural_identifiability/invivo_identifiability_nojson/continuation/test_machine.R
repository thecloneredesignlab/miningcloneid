# Labeled analytic fixtures only. No tumor model is evaluated.
test_machine<-function() {
 p<-c(log10_o2_S0=0);lo<-p;hi<-p+100
 start<-machine_new(p,0,lo,hi);trace<-list()
 f<-function(p,kind,step,direction,sequence){trace[[length(trace)+1L]]<<-list(p=p,kind=kind,step=step,direction=direction,sequence=sequence);-p[[1]]}
 no_observation<-function(p,value)NULL
 full<-machine_run(start,f,no_observation,100000L)
 stopifnot(full$reason=='poll_stationary',max(full$state$visits)>25L,all(full$state$visits>25L),full$state$best_p[[1]]>99.99)
 expected<-trace;trace<-list();state<-start
 repeat {
  part<-machine_run(state,f,no_observation,7L);state<-part$state
  if(part$reason=='poll_stationary')break
  stopifnot(part$reason=='call_limit')
 }
 stopifnot(identical(trace,expected),identical(state,full$state))
 timed<-machine_run(start,f,no_observation,10000L,deadline=0)
 stopifnot(timed$reason=='time_limit',timed$window_calls==0L,timed$state$total_calls==0L,timed$state$phase=='poll_evaluate')
 p<-c(log10_o2_S0=.1,log10_sigma_burden=0);lo<-p*0;hi<-lo+1
 value<-function(p)(p[[1]]-.8)^2+(p[[2]]-.5*p[[1]])^2
 observe<-function(p,v){q<-p;q[[2]]<-.5*q[[1]];list(parameters=q,expected=value(q))}
 start<-machine_new(p,value(p),lo,hi);trace<-list()
 evaluate<-function(p,kind,step,direction,sequence){trace[[length(trace)+1L]]<<-list(p=p,kind=kind,step=step,direction=direction,sequence=sequence);value(p)}
 whole<-machine_run(start,evaluate,observe,100000L);expected<-trace;trace<-list();state<-start;pending_seen<-FALSE
 repeat {
  part<-machine_run(state,evaluate,observe,1L);state<-part$state
  pending_seen<-pending_seen||state$phase=='observation_evaluate'
  if(part$reason=='poll_stationary')break
 }
 stopifnot(pending_seen,identical(trace,expected),identical(state,whole$state),state$observations>0L)
 list(pass=TRUE,beyond_25_visits=TRUE,all_scales_revisited=TRUE,partial_poll_sequence_identity=TRUE,
  observation_boundary_identity=TRUE,deadline_reason=TRUE,analytic_calls=whole$state$total_calls,
  maximum_visits=max(full$state$visits),model_calls=0L)
}
