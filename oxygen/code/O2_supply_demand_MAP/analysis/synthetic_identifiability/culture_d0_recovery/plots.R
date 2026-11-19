# Presentation only: base PDF devices require no X11 server or model calls.
d0_plot_recovery <- function(p,loss,out,branch) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 active<-p[as.logical(p$estimated),];stopifnot(nrow(active)==14L,all(is.finite(active$scaled_error)))
 draw<-function(name,width,height,fn) {
  grDevices::pdf(file.path(out,name),width=width,height=height,useDingbats=FALSE)
  on.exit(grDevices::dev.off());fn()
 }
 draw('parameter_recovery.pdf',11.25,7.1875,function(){
  par(mar=c(10,5,4,1));barplot(abs(active$scaled_error),names.arg=active$parameter,las=2,col='#3d6da8',
   ylab='Absolute error / registered natural parameter range',
   main=paste('Conditional D0 parameter error; origin',branch))
  mtext('14 free mechanisms; six nuisance values known. Descriptive errors, not intervals.',side=3,line=.4,cex=.8)
 })
 values<-as.numeric(unlist(loss[c('max_growth_error','max_kary_tv','max_flow_tv')]))
 stopifnot(length(values)==3L,all(is.finite(values)),all(values>=0))
 draw('prediction_recovery.pdf',8.75,5.625,function(){
  par(mar=c(5,5,4,1));barplot(values,names.arg=c('Growth / true sigma','Karyotype TV','Flow TV'),
   col='#3d6da8',ylim=c(0,max(.01,values)*1.15),ylab='Maximum discrepancy',
   main=paste('Expected observations; origin',branch))
  abline(h=.01,col='#af4035',lty=2);legend('topright','Registered target: each <= 0.01',lty=2,col='#af4035',bty='n',cex=.85)
 })
 invisible(NULL)
}
