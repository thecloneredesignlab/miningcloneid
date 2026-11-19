# Local-only adapter: use this study's registered resolution, not legacy bands.
d0_pattern <- function(fn, initial, directions, max_evaluations=20000L,
                       deadline=Inf, mesh=c(.02,.005,.001,.0001),
                       tolerance=1e-7, clock=function()as.numeric(Sys.time())) {
 n<-0L;exhausted<-FALSE;sweeps<-list();history<-list()
 evaluate<-function(x) {
  if(n>=max_evaluations||clock()>=deadline){exhausted<<-TRUE;return(Inf)}
  n<<-n+1L;q<-fn(x);if(length(q)!=1L||!is.finite(q))Inf else q
 }
 x<-pmin(pmax(initial,0),1);best<-evaluate(x)
 if(!is.finite(best)&&!exhausted)stop('Infeasible local start')
 if(is.finite(best))for(h in mesh) {
  repeat {
   before<-best;candidate<-x;candidate_q<-best;maximum<-0;complete<-TRUE
   for(i in seq_len(nrow(directions))) {
    for(sign in c(-1,1)) {
     z<-pmin(pmax(x+sign*h*directions[i,],0),1);if(identical(z,x))next
     q<-evaluate(z)
     if(exhausted){complete<-FALSE;break}
     maximum<-max(maximum,best-q)
     if(q<candidate_q){candidate<-z;candidate_q<-q}
    }
    if(exhausted)break
   }
   improvement<-best-candidate_q
   if(candidate_q<best){x<-candidate;best<-candidate_q}
   sweeps[[length(sweeps)+1L]]<-data.frame(mesh=h,complete=complete,before=before,after=best,improvement=improvement,maximum_tested_improvement=maximum,evaluations=n)
   history[[length(history)+1L]]<-x
   if(exhausted||improvement<=tolerance)break
  }
  if(exhausted)break
 }
 audit<-if(length(sweeps))do.call(rbind,sweeps)else data.frame()
 stable<-length(sweeps)>0L&&!exhausted&&all(vapply(tail(mesh,2),function(h){z<-audit[audit$mesh==h,,drop=FALSE];nrow(z)>0L&&isTRUE(tail(z$complete,1))&&tail(z$maximum_tested_improvement,1)<=tolerance},logical(1)))
 list(par=x,value=best,evaluations=n,sweeps=audit,tested_directional_stability=stable,budget_limited=exhausted,history=history,tolerance=tolerance)
}
d0_pattern_fixtures <- function() {
 q<-d0_pattern(function(x)sum(x*x),.004,matrix(1,1,1),100,mesh=c(.001,.0001))
 stopifnot(q$value<1e-16,q$tested_directional_stability,q$evaluations>4L)
 n<-0L;r<-d0_pattern(function(x){n<<-n+1L;sum(x*x)},rep(.9,14),diag(14),3)
 stopifnot(r$budget_limited,n==3L,r$evaluations==3L)
 t<-d0_pattern(function(x)stop('Must not evaluate after deadline'),.9,matrix(1,1,1),3,deadline=-1,clock=function()0)
 stopifnot(t$budget_limited,t$evaluations==0L)
 data.frame(fixture=c('Resolve improvements below legacy1e-5 to registered1e-7','Exact callback ceiling','Expired deadline before evaluation'),pass=TRUE,biological_calls=0L)
}
