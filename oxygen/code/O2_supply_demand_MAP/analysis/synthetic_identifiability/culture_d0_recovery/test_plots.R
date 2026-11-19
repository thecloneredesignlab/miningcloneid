# A headless fixture for the two presentation devices, with zero model calls.
source(commandArgs(TRUE)[1],local=.GlobalEnv);out<-commandArgs(TRUE)[2]
p<-data.frame(estimated=rep(TRUE,14),scaled_error=seq(0,.13,length.out=14),parameter=paste0('mechanism_',1:14))
d0_plot_recovery(p,list(max_growth_error=.001,max_kary_tv=.002,max_flow_tv=.003),out,'fixture')
for(n in c('parameter_recovery.pdf','prediction_recovery.pdf')) {
 f<-file.path(out,n);stopifnot(file.info(f)$size>1000,rawToChar(readBin(f,'raw',5))=='%PDF-')
}
cat('Headless PDF plotting fixture PASS; zero model calls\n')
