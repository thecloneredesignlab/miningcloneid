# Pure arithmetic oracles; no model calls. Also invoked inside the pinned SIF.
test_repair_arithmetic <- function() {
 set.seed(51931)
 # Enumerated stationary points must match an independent fine-grid oracle,
 # including active sigma limits, density bounds, volume floors and priors.
 for(K in c(0,.125))for(kind in c('ordinary','floor','all_floor')) {
  d<-list(y=c(-.4,.2,.7,-.1),logs=c(10,11,10.5,9.8),w=c(.125,.125,.375,.375),eps=1e-4,
   zlower=2,zupper=5,slo=.1,shi=.8,K=K,zmu=1.8)
  if(kind=='floor')d$logs<-c(-6,-5,-4,-7)
  if(kind=='all_floor')d$logs<-rep(-30,4)
  fit<-minimize_observation(d,3)
  grid<-seq(2,5,length.out=20001);values<-vapply(grid,function(z)observation_value(d,z),numeric(1))
  stopifnot(fit$value<=min(values)+1e-8,fit$z>=2,fit$z<=5,fit$sigma>=d$slo,fit$sigma<=d$shi)
  fixed<-minimize_observation(d,3,zfixed=4);stopifnot(fixed$z==4,abs(fixed$value-observation_value(d,4))<1e-12)
 }
 # Coupled nonsmooth valley has a jointly improving fixed block direction
 # where coordinate changes alone are unhelpful at the same poll scale.
 names_active<-c('log10_o2_S0','log10_kappa_O','log10_eta_o2')
 dirs<-pattern_directions(names_active);p<-setNames(c(.4,.4,.4),names_active)
 f<-function(q)abs(q[1]-q[2])+abs(q[2]-q[3])+.1*sum((q-.6)^2)
 coordinate<-vapply(dirs[seq_len(6)],function(v)f(p+.02*v),numeric(1))
 block<-vapply(dirs,function(v)f(p+.02*v),numeric(1))
 stopifnot(all(coordinate>=f(p)),min(block)<f(p))
 stopifnot(all(vapply(dirs,function(v)identical(names(v),names_active),logical(1))))
 # Distinct roots cannot recover through a shared winner assignment.
 roots<-c('archive_h25_production','formal_production_u1')
 stopifnot(length(unique(roots))==2L)
 TRUE
}
