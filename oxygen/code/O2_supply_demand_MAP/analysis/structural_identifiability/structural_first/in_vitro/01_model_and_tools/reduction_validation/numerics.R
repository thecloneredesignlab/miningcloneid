# Pure numerical audit helpers. No optimization or observation generation.
rr_counters<-c('native_buffer','biological','pruning','lower_out','upper_out','roundoff','error_events','ordinary_divisions','all_divisions','wgd','retained_loss','retained_gain','all_daughters')
rr_time<-seq(0,14,by=.25)
rr_rhs<-function(y,op){
 n<-length(op$grid);x<-y[seq_len(n)];L<-sum(x);crowd<-max(0,1-L/op$capacity)
 live<-crowd*as.numeric(op$G%*%x)-op$mu*x
 flow<-crowd*as.numeric(op$flows%*%x);direct<-sum(op$mu*x)
 q<-y[n+seq_len(7)];stock<-c(direct,flow[1:6])-op$clearance*q
 # Seven physical/audit stocks: direct death, native buffer, biological loss,
 # pruning, lower overflow, upper overflow, and signed roundoff residual.
 integrated<-c(direct,flow)
 c(live,stock,integrated)
}
rr_initial<-function(n0){c(as.numeric(n0),rep(0,7+1+length(rr_counters)))}
rr_euler<-function(op,n0,dt,times=rr_time){
 y<-rr_initial(n0);out<-matrix(0,length(times),length(y));out[1,]<-y;obs<-as.integer(round(times/dt));stopifnot(max(abs(obs*dt-times))<1e-10)
 clipped<-0;n<-length(op$grid);cursor<-2L
 for(k in seq_len(max(obs))){proposal<-y+dt*rr_rhs(y,op);negative<-proposal[seq_len(n)]<0;clipped<-clipped-sum(proposal[seq_len(n)][negative]);proposal[seq_len(n)]<-pmax(0,proposal[seq_len(n)]);y<-proposal
  if(cursor<=length(obs)&&k==obs[cursor]){out[cursor,]<-y;cursor<-cursor+1L}}
 list(times=times,state=out,op=op,method='Euler',dt=dt,clipped_live_equivalents=clipped)
}
rr_ode<-function(op,n0,rtol=1e-10,atol=1e-12,method='lsoda',times=rr_time){
 z<-deSolve::ode(y=rr_initial(n0),times=times,func=function(t,y,parameters)list(rr_rhs(y,op)),parms=NULL,method=method,rtol=rtol,atol=atol)
 state<-as.matrix(z[,-1,drop=FALSE]);n<-length(op$grid);negative<-min(state[,seq_len(n),drop=FALSE]);scale<-max(1,sum(n0));stopifnot(all(is.finite(state)),negative>=-1e-10*scale)
 list(times=times,state=state,op=op,method=method,rtol=rtol,atol=atol,minimum_live_state=negative)
}
rr_projection<-function(grid,bounds=c(32,54,76,98)){
 index<-vapply(grid,function(n)1L+sum(n>bounds),integer(1));outer(seq_len(length(bounds)+1L),index,'==')*1
}
rr_partition<-function(grid,candidate){
 stopifnot(min(grid)==22,max(grid)==154)
 original<-list(22:32,33:54,55:76,77:98,99:154);split_count<-switch(candidate,C5_original=1L,C5_uniform=1L,C10_uniform=2L,C20_uniform=4L,stop('Unknown candidate'))
 bins<-unlist(lapply(original,function(values){groups<-pmin(split_count,1L+floor((seq_along(values)-1L)*split_count/length(values)));split(values,groups)}),recursive=FALSE)
 P<-matrix(0,length(bins),length(grid));R<-matrix(0,length(grid),length(bins))
 for(j in seq_along(bins)){indices<-match(bins[[j]],grid);P[j,indices]<-1
  if(candidate=='C5_original')R[match(c(22,44,66,88,110)[j],grid),j]<-1 else R[indices,j]<-1/length(indices)}
 stopifnot(min(R)>=0,max(abs(P%*%R-diag(length(bins))))<1e-12)
 list(P=P,R=R,bins=bins,candidate=candidate)
}
rr_aggregate<-function(op,partition){
 P<-partition$P;R<-partition$R
 list(grid=as.numeric(colSums(op$grid*R)),G=as.matrix(P%*%op$G%*%R),mu=as.numeric(crossprod(R,op$mu)),flows=op$flows%*%R,capacity=op$capacity,clearance=op$clearance,lift=R,physical_grid=op$grid,candidate=partition$candidate)
}
rr_distribution<-function(run){
 n<-length(run$op$grid);live<-run$state[,seq_len(n),drop=FALSE]
 if(!is.null(run$op$lift)){live<-live%*%t(run$op$lift);grid<-run$op$physical_grid}else grid<-run$op$grid
 # Tiny solver negatives are visible in the raw state, but cannot define mass.
 live<-pmax(live,0);total<-rowSums(live);stopifnot(all(total>0))
 list(grid=grid,live=live,total=total,fraction=live/total,mean=as.numeric((live/total)%*%grid))
}
rr_compare<-function(a,b){
 stopifnot(identical(a$times,b$times));A<-rr_distribution(a);B<-rr_distribution(b);grid<-sort(unique(c(A$grid,B$grid)))
 pa<-matrix(0,nrow(A$fraction),length(grid));pb<-pa;pa[,match(A$grid,grid)]<-A$fraction;pb[,match(B$grid,grid)]<-B$fraction
 classP<-rr_projection(grid);class_tv<-rowSums(abs((pa-pb)%*%t(classP)))/2;full_tv<-rowSums(abs(pa-pb))/2
 arates<-rr_rates(a);brates<-rr_rates(b);mandatory<-c('direct','biological','error_events','wgd','retained_loss','retained_gain')
 gaps<-abs(arates[,mandatory,drop=FALSE]-brates[,mandatory,drop=FALSE]);limits<-.01*abs(arates[,mandatory,drop=FALSE])+1e-6
 c(max_class_tv=max(class_tv),max_full_tv=max(full_tv),max_log_live=max(abs(log(A$total/B$total))),max_mean_chromosome_gap=max(abs(A$mean-B$mean)),max_flux_tolerance_ratio=max(gaps/limits))
}
rr_rates<-function(run){
 n<-length(run$op$grid);x<-run$state[,seq_len(n),drop=FALSE];total<-rowSums(x);crowd<-pmax(0,1-total/run$op$capacity)
 flow<-t(run$op$flows%*%t(x));flow<-flow*crowd/total;colnames(flow)<-rr_counters
 cbind(direct=as.numeric(x%*%run$op$mu)/total,flow)
}
rr_gate_state<-function(v){v[['max_class_tv']]<=.05&&v[['max_log_live']]<=.1&&v[['max_mean_chromosome_gap']]<=1&&v[['max_flux_tolerance_ratio']]<=1}
rr_gate_numeric<-function(v,tolerance=1e-6,mean_tolerance=1e-4){v[['max_full_tv']]<=tolerance&&v[['max_log_live']]<=tolerance&&v[['max_mean_chromosome_gap']]<=mean_tolerance}
rr_moments<-function(run){
 d<-rr_distribution(run);P<-rr_projection(d$grid);rows<-list()
 for(j in seq_len(nrow(P))){idx<-which(P[j,]>0);counts<-rowSums(d$live[,idx,drop=FALSE]);mean<-as.numeric(d$live[,idx,drop=FALSE]%*%d$grid[idx])/counts;second<-as.numeric(d$live[,idx,drop=FALSE]%*%d$grid[idx]^2)/counts
  mean[counts==0]<-NA_real_;second[counts==0]<-NA_real_;rows[[j]]<-data.frame(time=run$times,class=j,live=counts,mean_chromosomes=mean,variance_chromosomes=pmax(0,second-mean^2))}
 do.call(rbind,rows)
}
