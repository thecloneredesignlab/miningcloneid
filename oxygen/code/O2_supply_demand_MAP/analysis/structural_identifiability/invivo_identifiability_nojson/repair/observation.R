# Conditional observation-block proposals; the canonical scorer accepts every move.
observation_data <- function(ctx,record) {
 y<-logs<-w<-numeric();H<-length(ctx$scenarios)
 stopifnot(H==8L,length(record$predictions$native)==H)
 for(i in seq_len(H)) {
  sc<-ctx$scenarios[[i]];n<-record$predictions$native[[i]]$result
  keep<-ctx$cfg$scenario_cpp$keep_burden[[i]] & is.finite(sc$obs_burden) & sc$obs_burden>0
  stock<-n$Ntot_total_obs[keep];obs<-sc$obs_burden[keep]
  stopifnot(length(stock)==length(obs),length(obs)>0L,all(is.finite(stock)),all(stock>0))
  y<-c(y,log(pmax(obs,ctx$cfg$burden_log_eps)));logs<-c(logs,log(stock));w<-c(w,rep(1/(H*length(obs)),length(obs)))
 }
 stopifnot(length(y)==92L,abs(sum(w)-1)<1e-14)
 list(y=y,logs=logs,w=w,eps=ctx$cfg$burden_log_eps,
  zlower=log(10)*ctx$lower[['log10_rho_2N']],zupper=log(10)*ctx$upper[['log10_rho_2N']],
  slo=10^ctx$lower[['log10_sigma_burden']],shi=10^ctx$upper[['log10_sigma_burden']],
  K=if(ctx$policy=="production") ctx$cfg$lambda_prior/(2*ctx$cfg$prior_sd_log10_rho_2N^2*log(10)^2) else 0,
  zmu=log(10)*ctx$cfg$prior_center_log10_rho_2N)
}
observation_A <- function(d,z) sum(d$w*(d$y-pmax(d$logs-z,log(d$eps)))^2)
observation_value <- function(d,z,sigma=NULL) {
 A<-observation_A(d,z)
 if(is.null(sigma)) sigma<-min(d$shi,max(d$slo,sqrt(max(0,A))))
 log(sigma)+.5*log(2*pi)+A/(2*sigma^2)+d$K*(z-d$zmu)^2
}
real_polynomial_roots <- function(coeff) {
 while(length(coeff)>1L && abs(tail(coeff,1))<1e-16) coeff<-head(coeff,-1L)
 if(length(coeff)<2L) return(numeric())
 z<-polyroot(coeff);Re(z[abs(Im(z))<1e-7*(1+abs(Re(z)))])
}
minimize_observation <- function(d,zcurrent,zfixed=NULL) {
 stopifnot(d$zlower<=d$zupper,d$slo>0,d$shi>=d$slo)
 candidates<-c(d$zlower,d$zupper,zcurrent)
 if(!is.null(zfixed)) candidates<-zfixed else {
  breaks<-sort(unique(c(d$zlower,d$zupper,(d$logs-log(d$eps))[d$logs-log(d$eps)>d$zlower & d$logs-log(d$eps)<d$zupper])))
  candidates<-c(candidates,breaks)
  for(i in seq_len(length(breaks)-1L)) {
   lo<-breaks[i];hi<-breaks[i+1L];mid<-(lo+hi)/2;active<-d$logs-mid>log(d$eps)
   delta<-d$y-d$logs;a<-sum(d$w[active]);b<-2*sum(d$w[active]*delta[active])
   cc<-sum(d$w[active]*delta[active]^2)+sum(d$w[!active]*(d$y[!active]-log(d$eps))^2)
   roots<-numeric()
   for(sig in c(d$slo,d$shi)) {
    den<-a/sig^2+2*d$K
    if(den>0) roots<-c(roots,(2*d$K*d$zmu-b/(2*sig^2))/den)
    roots<-c(roots,real_polynomial_roots(c(cc-sig^2,b,a)))
   }
   roots<-c(roots,real_polynomial_roots(c(b-4*d$K*d$zmu*cc,2*a+4*d$K*(cc-d$zmu*b),4*d$K*(b-d$zmu*a),4*d$K*a)))
   candidates<-c(candidates,roots[is.finite(roots) & roots>=lo & roots<=hi])
  }
 }
 candidates<-sort(unique(candidates[is.finite(candidates) & candidates>=d$zlower-1e-12 & candidates<=d$zupper+1e-12]))
 stopifnot(length(candidates)>0L)
 values<-vapply(candidates,function(z) observation_value(d,z),numeric(1));best<-min(values)
 ties<-which(values<=best+1e-13);idx<-ties[which.min(abs(candidates[ties]-zcurrent))]
 z<-candidates[idx];sigma<-min(d$shi,max(d$slo,sqrt(max(0,observation_A(d,z)))))
 list(z=z,sigma=sigma,value=observation_value(d,z,sigma),candidates=candidates,values=values)
}
observation_proposal <- function(ctx,record,fixed=NULL) {
 d<-observation_data(ctx,record);p<-record$parameters
 oldz<-log(10)*p[['log10_rho_2N']];olds<-10^p[['log10_sigma_burden']]
 oldburden<-observation_value(d,oldz,olds)-d$K*(oldz-d$zmu)^2
 stopifnot(abs(oldburden-record$components$L_b)<1e-8)
 fit<-minimize_observation(d,oldz,if('log10_rho_2N'%in%names(fixed)) log(10)*fixed[['log10_rho_2N']] else NULL)
 q<-p
 q[['log10_rho_2N']]<-if('log10_rho_2N'%in%names(fixed)) fixed[['log10_rho_2N']] else min(ctx$upper[['log10_rho_2N']],max(ctx$lower[['log10_rho_2N']],fit$z/log(10)))
 q[['log10_sigma_burden']]<-min(ctx$upper[['log10_sigma_burden']],max(ctx$lower[['log10_sigma_burden']],log10(fit$sigma)))
 z<-log(10)*q[['log10_rho_2N']];sigma<-10^q[['log10_sigma_burden']]
 expected<-record$components$L+observation_value(d,z,sigma)-observation_value(d,oldz,olds)
 stopifnot(expected<=record$components$L+1e-8,all(q[names(fixed)]==fixed),all(q>=ctx$lower),all(q<=ctx$upper))
 list(parameters=q,expected=expected,expected_burden=observation_value(d,z,sigma)-d$K*(z-d$zmu)^2,
  candidate_count=length(fit$candidates),source_objective=record$components$L)
}
