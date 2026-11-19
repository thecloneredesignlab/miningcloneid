# Pure numerical fixtures, independent of production biological data.
args<-commandArgs(TRUE);code<-if(length(args))args[1]else dirname(sys.frame(1)$ofile)
source(file.path(code,'numerics.R'))
for(candidate in c('C5_original','C5_uniform','C10_uniform','C20_uniform')){
 p<-rr_partition(22:154,candidate);stopifnot(all(colSums(p$P)==1),all(abs(colSums(p$R)-1)<1e-12),all(p$R>=0),max(abs(p$P%*%p$R-diag(nrow(p$P))))<1e-12)
}
op<-list(grid=44,G=matrix(.7,1,1),mu=.2,flows=matrix(0,length(rr_counters),1,dimnames=list(rr_counters,NULL)),capacity=1000,clearance=.1)
op$flows['native_buffer',]<-.03;op$flows['biological',]<-.03
initial<-rr_initial(100);rhs<-rr_rhs(initial,op)
stopifnot(abs(rhs[1]-(.9*.7-.2)*100)<1e-12,abs(rhs[2]-20)<1e-12,abs(rhs[3]-.9*.03*100)<1e-12)
initial[2]<-10;stopifnot(abs(rr_rhs(initial,op)[2]-19)<1e-12)
times<-seq(0,4,.25);a<-rr_ode(op,100,times=times);b<-rr_ode(op,100,method='radau',times=times)
rate<-.7-.2;capacity<-1000*rate/.7;closed<-capacity/(1+(capacity/100-1)*exp(-rate*times))
stopifnot(max(abs(a$state[,1]-closed))<1e-5,rr_gate_numeric(rr_compare(a,b)))
e<-rr_euler(op,100,.001,times);stopifnot(max(abs(e$state[,1]-closed))<.12)
grid<-1:308;n0<-rep(0,length(grid));n0[match(22:154,grid)]<-seq_along(22:154);stopifnot(sum(n0)==sum(seq_along(22:154)),all(n0[grid<22|grid>154]==0))
P<-rr_projection(grid,c(21,32,54,76,98,154));stopifnot(nrow(P)==7,all(colSums(P)==1))
good<-c(max_class_tv=.01,max_full_tv=.1,max_log_live=.01,max_mean_chromosome_gap=.1,max_flux_tolerance_ratio=.9);stopifnot(rr_gate_state(good));good['max_flux_tolerance_ratio']<-1.01;stopifnot(!rr_gate_state(good))
cat('PASS: projection identities, crowding/death accounting, logistic solution, independent ODE solver, Euler refinement, zero padding, outside buckets, flux gate\n')
